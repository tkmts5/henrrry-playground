/**
 * rough-three.js
 * Three.js 用のラフ描画ライブラリ（数学座標系 x=横, y=奥行き, z=高さ）。
 * THREE はグローバルに読み込まれた状態で使用すること。
 *
 * 使用例:
 *   <script src="three.min.js"></script>
 *   <script src="rough-three.js"></script>
 *   RoughThree.drawRoughFilledBox(scene, camera, center, width, depth, height);
 *
 * drawRoughFilledSphere: ラフな球体を描画
 */
(function (global) {
    'use strict';

    const THREE = global.THREE;
    if (!THREE) {
        console.error('rough-three.js: THREE が読み込まれていません。');
        return;
    }

    const DEFAULT_BOX_COLOR = '#f5edd8';
    const DEFAULT_SPHERE_COLOR = '#f5edd8';

    const PERLIN_P = new Uint8Array(512);
    (function () {
        const ref = [151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20,125,136,171,168,68,175,74,165,71,134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,102,143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200,196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226,250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,223,183,170,213,119,248,152,2,44,154,163,70,221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,224,232,178,185,112,104,218,246,97,228,251,31,177,127,193,211,230,241,235,249,14,239,107,49,192,214,31,181,199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180];
        for (let i = 0; i < 256; i++) PERLIN_P[i] = ref[i];
        for (let i = 256; i < 512; i++) PERLIN_P[i] = PERLIN_P[i & 255];
    })();
    const PERLIN_G3 = [
        [1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0],
        [1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1],
        [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1]
    ];
    function perlinFade(t) {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }
    function perlinLerp(a, b, t) {
        return a + t * (b - a);
    }
    //const PERLIN_EPS = 1e-10;
    function perlinNoise3D(x, y, z) {
        const X = Math.floor(x) & 255, Y = Math.floor(y) & 255, Z = Math.floor(z) & 255;
        x -= Math.floor(x);
        y -= Math.floor(y);
        z -= Math.floor(z);
        //if (x === 0) x = PERLIN_EPS;
        //if (y === 0) y = PERLIN_EPS;
        //if (z === 0) z = PERLIN_EPS;
        const u = perlinFade(x), v = perlinFade(y), w = perlinFade(z);
        const A = (PERLIN_P[X] + Y) & 255, B = (PERLIN_P[X + 1] + Y) & 255;
        const AA = (PERLIN_P[A] + Z) & 255, AB = (PERLIN_P[A + 1] + Z) & 255;
        const BA = (PERLIN_P[B] + Z) & 255, BB = (PERLIN_P[B + 1] + Z) & 255;
        const grad = (h, dx, dy, dz) => {
            const g = PERLIN_G3[h % 12];
            return g[0] * dx + g[1] * dy + g[2] * dz;
        };
        const AA1 = (AA + 1) & 255, BA1 = (BA + 1) & 255, AB1 = (AB + 1) & 255, BB1 = (BB + 1) & 255;
        return perlinLerp(
            perlinLerp(
                perlinLerp(grad(PERLIN_P[AA], x, y, z), grad(PERLIN_P[BA], x - 1, y, z), u),
                perlinLerp(grad(PERLIN_P[AB], x, y - 1, z), grad(PERLIN_P[BB], x - 1, y - 1, z), u),
                v
            ),
            perlinLerp(
                perlinLerp(grad(PERLIN_P[AA1], x, y, z - 1), grad(PERLIN_P[BA1], x - 1, y, z - 1), u),
                perlinLerp(grad(PERLIN_P[AB1], x, y - 1, z - 1), grad(PERLIN_P[BB1], x - 1, y - 1, z - 1), u),
                v
            ),
            w
        );
    }
    function perlinNoise2D(x, y) {
        return perlinNoise3D(x, y, 0);
    }

    /**
     * 数学座標 (x=横, y=奥行き, z=高さ) を Three.js 座標に変換する。
     */
    function mathToThree(mx, my, mz) {
        return new THREE.Vector3(mx, mz, -my);
    }

    /**
     * 辺 a→b を segmentLength ごとに区切り、中間点にのみ jitter を加えてラフな点列を返す。
     * 端点は a, b をそのまま使い、3辺が交わる頂点でずれが生じないようにする。
     * edgeSeed を渡すと Perlin ノイズで決定的なゆらぎになり、同じ edgeSeed なら常に同じパターンになる
     * （ボクセルアートで隣接ボックスの隙間を防ぐため）。
     * @param {THREE.Vector3} a - 始点
     * @param {THREE.Vector3} b - 終点
     * @param {number} segmentLength - 区切る長さ（約何pxごと）
     * @param {number} jitter - ずれの大きさ
     * @param {number} [edgeSeed] - 0,1,2 など。指定時は決定的なラフ（未指定時は Math.random）
     * @returns {THREE.Vector3[]}
     */
    function roughEdgePoints(a, b, segmentLength, jitter, edgeSeed) {
        const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
        const totalLen = Math.hypot(dx, dy, dz);
        const steps = Math.max(1, Math.ceil(totalLen / segmentLength));
        const points = [];
        const useDeterministic = typeof edgeSeed === 'number';
        for (let i = 0; i <= steps; i++) {
            const t = i / steps;
            if (i === 0) {
                points.push(a.clone());
            } else if (i === steps) {
                points.push(b.clone());
            } else {
                let jx, jy, jz;
                if (useDeterministic) {
                    const u = i * 1.23, v = edgeSeed * 71;
                    jx = perlinNoise3D(u, v, 0) * jitter;
                    jy = perlinNoise3D(u + 50, v, 0.3) * jitter;
                    jz = perlinNoise3D(u + 100, v, 0.7) * jitter;
                } else {
                    jx = (Math.random() - 0.5) * 2 * jitter;
                    jy = (Math.random() - 0.5) * 2 * jitter;
                    jz = (Math.random() - 0.5) * 2 * jitter;
                }
                points.push(new THREE.Vector3(
                    a.x + dx * t + jx,
                    a.y + dy * t + jy,
                    a.z + dz * t + jz
                ));
            }
        }
        return points;
    }

    const BOX_EDGES = [
        [0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7],
        [4, 5], [4, 6], [5, 7], [6, 7]
    ];

    /** 向かい合う4本ずつで同じ形状にするための辺の分類（幅・奥行き・高さの3種類） */
    const BOX_EDGE_TYPES = {
        width:  [[0, 4], [1, 5], [2, 6], [3, 7]],
        depth:  [[0, 2], [1, 3], [4, 6], [5, 7]],
        height: [[0, 1], [2, 3], [4, 5], [6, 7]]
    };

    /**
     * 1本の辺から (t, jitter) のテンプレートを生成する。同じ形状を平行な辺に適用するため。
     * @param {THREE.Vector3[]} points - roughEdgePoints の戻り値
     * @param {THREE.Vector3} a - 始点
     * @param {THREE.Vector3} b - 終点
     * @returns {{ t: number, jx: number, jy: number, jz: number }[]}
     */
    function edgeTemplateFromPoints(points, a, b) {
        const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
        const lenSq = dx * dx + dy * dy + dz * dz || 1;
        return points.map((p) => {
            const px = p.x - a.x, py = p.y - a.y, pz = p.z - a.z;
            const t = (px * dx + py * dy + pz * dz) / lenSq;
            return {
                t,
                jx: p.x - (a.x + dx * t),
                jy: p.y - (a.y + dy * t),
                jz: p.z - (a.z + dz * t)
            };
        });
    }

    /**
     * テンプレートを別の辺 A'→B' に適用し、同じぶらし形状の点列を返す。
     */
    function applyEdgeTemplate(vertices, lo, hi, template) {
        const a = vertices[lo], b = vertices[hi];
        const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
        return template.map(({ t, jx, jy, jz }) => new THREE.Vector3(
            a.x + dx * t + jx,
            a.y + dy * t + jy,
            a.z + dz * t + jz
        ));
    }

    /**
     * 数学座標系で直方体の6面をラフな輪郭で生成し、sceneに追加する。
     * 8頂点を共有し、12本の辺を1回だけ計算して接する面で共有する。
     * カメラから見て裏を向いている面は描画しない。
     *
     * @param {THREE.Scene} scene
     * @param {THREE.Camera} camera - 可視判定に使用
     * @param {THREE.Vector3} center - 中心（数学座標 x,y,z）
     * @param {number} width  - 横（x軸方向）
     * @param {number} depth  - 縦・奥行き（y軸方向）
     * @param {number} height - 高さ（z軸方向）
     * @param {string} [color] - ベース色（省略時は DEFAULT_BOX_COLOR）
     * @param {number} [segmentLength] - 辺を区切る長さ
     * @param {number} [jitter] - ずれの大きさ
     */
    function drawRoughFilledBox(scene, camera, center, width, depth, height, color, segmentLength, jitter) {
        const segLen = segmentLength ?? 3;//3
        const jit = jitter ?? 0.2;//0.2
        const baseColor = color ?? DEFAULT_BOX_COLOR;
        const wx = width / 2, dy = depth / 2, hz = height / 2;
        const mathVertices = [
            [center.x + wx, center.y + dy, center.z + hz],
            [center.x + wx, center.y + dy, center.z - hz],
            [center.x + wx, center.y - dy, center.z + hz],
            [center.x + wx, center.y - dy, center.z - hz],
            [center.x - wx, center.y + dy, center.z + hz],
            [center.x - wx, center.y + dy, center.z - hz],
            [center.x - wx, center.y - dy, center.z + hz],
            [center.x - wx, center.y - dy, center.z - hz]
        ];
        const vertices = mathVertices.map(([mx, my, mz]) => mathToThree(mx, my, mz));

        // 3種類の辺形状テンプレートを、各代表辺1本だけ rough で生成
        const [canonX, canonY, canonZ] = [
            BOX_EDGE_TYPES.width[0],
            BOX_EDGE_TYPES.depth[0],
            BOX_EDGE_TYPES.height[0]
        ];
        const templateX = edgeTemplateFromPoints(
            roughEdgePoints(vertices[canonX[0]], vertices[canonX[1]], segLen, jit, 0),
            vertices[canonX[0]], vertices[canonX[1]]
        );
        const templateY = edgeTemplateFromPoints(
            roughEdgePoints(vertices[canonY[0]], vertices[canonY[1]], segLen, jit, 1),
            vertices[canonY[0]], vertices[canonY[1]]
        );
        const templateZ = edgeTemplateFromPoints(
            roughEdgePoints(vertices[canonZ[0]], vertices[canonZ[1]], segLen, jit, 2),
            vertices[canonZ[0]], vertices[canonZ[1]]
        );

        const edgeKeyToTemplate = new Map();
        BOX_EDGE_TYPES.width.forEach(([p, q]) => edgeKeyToTemplate.set(`${Math.min(p, q)}-${Math.max(p, q)}`, templateX));
        BOX_EDGE_TYPES.depth.forEach(([p, q]) => edgeKeyToTemplate.set(`${Math.min(p, q)}-${Math.max(p, q)}`, templateY));
        BOX_EDGE_TYPES.height.forEach(([p, q]) => edgeKeyToTemplate.set(`${Math.min(p, q)}-${Math.max(p, q)}`, templateZ));

        const edgeCache = new Map();
        BOX_EDGES.forEach(([i, j]) => {
            const lo = Math.min(i, j), hi = Math.max(i, j);
            const template = edgeKeyToTemplate.get(`${lo}-${hi}`);
            edgeCache.set(`${lo}-${hi}`, applyEdgeTemplate(vertices, lo, hi, template));
        });

        function getEdgePoints(fromIdx, toIdx) {
            const lo = Math.min(fromIdx, toIdx), hi = Math.max(fromIdx, toIdx);
            const points = edgeCache.get(`${lo}-${hi}`);
            return fromIdx < toIdx ? points : [...points].reverse();
        }

        const faceList = [
            [0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]
        ];

        const material = new THREE.MeshLambertMaterial({
            color: baseColor,
            side: THREE.DoubleSide
        });

        const cameraPos = camera.position.clone();

        faceList.forEach((indices) => {
            const [a, b, c, d] = indices;
            const va = vertices[a], vb = vertices[b], vc = vertices[c], vd = vertices[d];
            const faceCenter = new THREE.Vector3().addVectors(va, vb).add(vc).add(vd).divideScalar(4);
            const faceNormal = new THREE.Vector3().crossVectors(
                new THREE.Vector3().subVectors(vb, va),
                new THREE.Vector3().subVectors(vd, va)
            ).normalize();
            const toCamera = new THREE.Vector3().subVectors(cameraPos, faceCenter);
            if (faceNormal.dot(toCamera) <= 0) return;

            const e0 = getEdgePoints(a, b);
            const e1 = getEdgePoints(b, c);
            const e2 = getEdgePoints(c, d);
            const e3 = getEdgePoints(d, a);
            // 8頂点は座標をぶらさず、vertices の座標をそのまま使う（e0/e1/e2/e3 の端点は中間点のみ使う）
            const outline = [
                vertices[a].clone(),
                ...e0.slice(1, -1),
                vertices[b].clone(),
                ...e1.slice(1, -1),
                vertices[c].clone(),
                ...e2.slice(1, -1),
                vertices[d].clone(),
                ...e3.slice(1, -1)
            ];
            const n = outline.length;
            const centroid = new THREE.Vector3(0, 0, 0);
            outline.forEach((p) => centroid.add(p));
            centroid.divideScalar(n);
            const positions = [centroid.x, centroid.y, centroid.z];
            outline.forEach((p) => { positions.push(p.x, p.y, p.z); });
            const geometry = new THREE.BufferGeometry();
            geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(positions), 3));
            const triangleIndices = [];
            for (let i = 0; i < n; i++) {
                triangleIndices.push(0, 1 + i, 1 + ((i + 1) % n));
            }
            geometry.setIndex(triangleIndices);
            geometry.computeVertexNormals();
            scene.add(new THREE.Mesh(geometry, material));
        });
    }

    function icosphereVerticesAndFaces(radius, detail) {
        const t = (1 + Math.sqrt(5)) / 2;
        const verts = [
            [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
            [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
            [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]
        ].map(([x, y, z]) => {
            const r = Math.hypot(x, y, z);
            return new THREE.Vector3(x / r * radius, y / r * radius, z / r * radius);
        });
        let faces = [
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]
        ];
        const edgeKey = (a, b) => (a < b ? a + ',' + b : b + ',' + a);
        for (let d = 0; d < detail; d++) {
            const mid = new Map();
            const nextFaces = [];
            faces.forEach(([a, b, c]) => {
                let mabI = mid.get(edgeKey(a, b));
                if (mabI === undefined) {
                    mabI = verts.length;
                    verts.push(verts[a].clone().add(verts[b]).multiplyScalar(0.5));
                    mid.set(edgeKey(a, b), mabI);
                }
                let mbcI = mid.get(edgeKey(b, c));
                if (mbcI === undefined) {
                    mbcI = verts.length;
                    verts.push(verts[b].clone().add(verts[c]).multiplyScalar(0.5));
                    mid.set(edgeKey(b, c), mbcI);
                }
                let mcaI = mid.get(edgeKey(c, a));
                if (mcaI === undefined) {
                    mcaI = verts.length;
                    verts.push(verts[c].clone().add(verts[a]).multiplyScalar(0.5));
                    mid.set(edgeKey(c, a), mcaI);
                }
                nextFaces.push([a, mabI, mcaI], [b, mbcI, mabI], [c, mcaI, mbcI], [mabI, mbcI, mcaI]);
            });
            faces = nextFaces;
            for (let i = 12; i < verts.length; i++) {
                verts[i].normalize().multiplyScalar(radius);
            }
        }
        return { verts, faces };
    }

    /**
     * 数学座標系でラフな球体を生成し、sceneに追加する。
     * ICO球を頂点共有で自前構築し、全頂点の座標を定めてから jitter を1回ずつ適用する。
     */
    function drawRoughFilledSphere(scene, center, radius, color, segmentLength, jitter) {
        const segLen = segmentLength ?? 2;//Math.max(4, radius * 0.3);
        const jit = jitter ?? 0.1;//radius * 0.03;
        const baseColor = color ?? DEFAULT_SPHERE_COLOR;
        const edgeApprox = radius * 1.05;
        const detail = Math.max(0, Math.min(4, Math.ceil(Math.log2(edgeApprox / segLen))));
        const { verts, faces } = icosphereVerticesAndFaces(radius, detail);
        verts.forEach((v) => {
            const r = Math.hypot(v.x, v.y, v.z) || 1;
            const j = (Math.random() - 0.5) * 2 * jit;
            v.x += (v.x / r) * j;
            v.y += (v.y / r) * j;
            v.z += (v.z / r) * j;
        });
        const positions = new Float32Array(verts.length * 3);
        verts.forEach((v, i) => {
            positions[i * 3] = v.x;
            positions[i * 3 + 1] = v.y;
            positions[i * 3 + 2] = v.z;
        });
        const indices = faces.flat();
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geometry.setIndex(indices);
        geometry.computeVertexNormals();
        const material = new THREE.MeshLambertMaterial({
            color: baseColor,
            side: THREE.DoubleSide
        });
        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.copy(mathToThree(center.x, center.y, center.z));
        scene.add(mesh);
    }

    // -------------------------------------------------------------------------
    // ボクセルアート（7×7×7、129文字 + 256文字変換表 + パレット）
    // -------------------------------------------------------------------------
    const VOXEL_SIZE = 7;
    const VOXEL_COUNT = VOXEL_SIZE * VOXEL_SIZE * VOXEL_SIZE;
    const VOXEL_PACKED_LEN = Math.ceil((VOXEL_COUNT * 3) / 8);
    const DEFAULT_VOXEL_EXPORT_TABLE = "ぁあぃいぅうぇえぉおかがきぎくぐけげこごさざしじすずせぜそぞただちぢっつづてでとどなにぬねのはばぱひびぴふぶぷへべぺほぼぽまみむめもゃやゅゆょよらりるれろゎわゐゑをんゔゕゖァアィイゥウェエォオカガキギクグケゲコゴサザシジスズセゼソゾタダチヂッツヅテデトドナニヌネノハバパヒビピフブプヘベペホボポマミムメモャヤュユョヨラリルレロヮワヲンヴヵヶ一右雨円王音下火花貝学気九休玉金空月犬見五口校左三山子四糸字耳七車手十出女小上森人水正生青夕石赤千川先早草足村大男竹中虫町天田土二日入年白八百文木本名目立力林六引羽雲園遠何";

    /**
     * 129文字のデータ文字列を343要素のボクセル配列（各要素0〜7）に復元する。
     * 各文字は変換表のインデックス（0〜255）でバイト値に変換し、8ビットずつ並べて1029ビット、
     * さらに3ビットずつ343個の値に展開する。
     * 変換表は引数またはグローバル VOXEL_EXPORT_TABLE がなければ、組み込みの DEFAULT_VOXEL_EXPORT_TABLE を使用する。
     * @param {string} s - 129文字の文字列（voxel-edit のエクスポート形式）
     * @param {string} [exportTable] - 256文字の変換表。省略時はグローバルまたは組み込み表を使用
     * @returns {Uint8Array} 長さ343、各要素0〜7
     */
    function unpackVoxelData(s, exportTable) {
        const g = typeof global !== 'undefined' ? global : (typeof window !== 'undefined' ? window : {});
        const table = exportTable || g.VOXEL_EXPORT_TABLE || DEFAULT_VOXEL_EXPORT_TABLE;
        const voxels = new Uint8Array(VOXEL_COUNT);
        if (table.length < 256) return voxels;
        var str = s;
        if (str.length === VOXEL_PACKED_LEN - 1) str = str + table[0];
        if (str.length !== VOXEL_PACKED_LEN) return voxels;
        const bits = [];
        for (let j = 0; j < VOXEL_PACKED_LEN; j++) {
            const idx = table.indexOf(str[j]);
            const byteVal = (idx >= 0 ? idx : 0) & 255;
            for (let bitInByte = 0; bitInByte < 8; bitInByte++) bits.push((byteVal >> bitInByte) & 1);
        }
        for (let i = 0; i < VOXEL_COUNT; i++) {
            let v = 0;
            for (let b = 0; b < 3; b++) if (bits[i * 3 + b]) v |= 1 << b;
            voxels[i] = v;
        }
        return voxels;
    }

    /**
     * 129文字のデータ文字列とカラーパレットを受け取り、Three.js の scene 上に
     * 7×7×7 のボクセルアートを drawRoughFilledBox で描画する。
     * 数学座標系: x=横, y=奥行き, z=高さ。drawRoughFilledBox と同じく center は数学座標の THREE.Vector3。
     * 変換表は組み込み済み。voxel-export-table.js を読み込めばその表を優先使用。
     *
     * @param {THREE.Scene} scene
     * @param {THREE.Camera} camera - drawRoughFilledBox の可視判定に使用
     * @param {string} dataString - 129文字の文字列（voxel-edit エクスポート形式、256文字変換表対応）
     * @param {string[]} palette - 色の配列。palette[0] は未使用、palette[1]〜palette[7] が色（例: '#ff0000'）
     * @param {THREE.Vector3} [center] - ボクセルアート全体の中心（数学座標）。省略時は原点
     * @param {number} [boxSize=1] - 1ボックスの1辺のサイズ（数学座標）
     * @param {number} [segmentLength] - drawRoughFilledBox に渡す辺の区切り長さ
     * @param {number} [jitter] - drawRoughFilledBox に渡すずれの大きさ
     */
    function drawVoxelArt(scene, camera, dataString, palette, center, boxSize, segmentLength, jitter) {
        var len = typeof dataString === 'string' ? dataString.length : 0;
        if (len !== VOXEL_PACKED_LEN && len !== VOXEL_PACKED_LEN - 1) {
            console.warn('rough-three.js drawVoxelArt: dataString は' + VOXEL_PACKED_LEN + '文字（または' + (VOXEL_PACKED_LEN - 1) + '文字）である必要があります。現在: ' + len);
            return;
        }
        const hasCenter = center && center.isVector3;
        const artCenter = hasCenter ? center : new THREE.Vector3(0, 0, 0);
        const size = hasCenter
            ? ((boxSize != null && boxSize > 0) ? Number(boxSize) : 1)
            : ((center != null && Number(center) > 0) ? Number(center) : 1);
        const segLen = hasCenter ? segmentLength : boxSize;
        const jit = hasCenter ? jitter : segmentLength;
        const voxels = unpackVoxelData(dataString);
        const SIZE = VOXEL_SIZE;
        const half = (SIZE / 2) * size;
        for (let z = 0; z < SIZE; z++) {
            for (let y = 0; y < SIZE; y++) {
                for (let x = 0; x < SIZE; x++) {
                    const c = voxels[x + SIZE * y + SIZE * SIZE * z];
                    if (c === 0) continue;
                    const color = palette && palette[c] ? palette[c] : DEFAULT_BOX_COLOR;
                    const mx = (x + 0.5) * size - half + artCenter.x;
                    const my = ((SIZE - 1 - z) + 0.5) * size - half + artCenter.y;
                    const mz = (y + 0.5) * size - half + artCenter.z;
                    const boxCenter = new THREE.Vector3(mx, my, mz);
                    drawRoughFilledBox(scene, camera, boxCenter, size, size, size, color, segLen, jit);
                    /*drawRoughFilledSphere(
                        scene,
                        boxCenter,
                        size, color,
                        segLen, jit
                    );*/
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // 公開 API
    // -------------------------------------------------------------------------
    global.RoughThree = {
        mathToThree,
        roughEdgePoints,
        drawRoughFilledBox,
        drawRoughFilledSphere,
        drawVoxelArt,
        unpackVoxelData,
        perlinNoise2D,
        perlinNoise3D,
        DEFAULT_BOX_COLOR,
        DEFAULT_SPHERE_COLOR,
        VOXEL_PACKED_LEN,
        VOXEL_SIZE,
        VOXEL_COUNT
    };

})(typeof window !== 'undefined' ? window : this);
