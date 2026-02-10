// wormhole-viz.js - Three.js visualization with clean KaTeX derivations

// ============================================================================
// INITIALIZATION
// ============================================================================

const canvas = document.getElementById('wormhole-canvas');
const container = document.getElementById('canvas-container');

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0a0a0f);

const camera = new THREE.PerspectiveCamera(60, container.clientWidth / container.clientHeight, 0.1, 1000);
camera.position.set(8, 5, 8);
camera.lookAt(0, 0, 0);

const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setSize(container.clientWidth, container.clientHeight);
renderer.setPixelRatio(window.devicePixelRatio);

// Lighting
const ambientLight = new THREE.AmbientLight(0x404060, 0.5);
scene.add(ambientLight);

const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
directionalLight.position.set(5, 10, 5);
scene.add(directionalLight);

const pointLight = new THREE.PointLight(0x7799ff, 0.5, 20);
pointLight.position.set(0, 0, 0);
scene.add(pointLight);

// ============================================================================
// WORMHOLE PARAMETERS & GEOMETRY
// ============================================================================

let wormholeMesh = null;
let gridHelper = null;
let geodesicLines = [];

let params = {
    shapeType: 'constant',
    b0: 1.0,
    vizMode: 'embedding'
};

function getShapeConfig() {
    return TensorCalc.WormholeConfigs[params.shapeType](params.b0);
}

function embeddingZ(r, config) {
    const b0 = params.b0;
    if (r <= b0) return 0;
    
    const steps = 100;
    const dr = (r - b0 * 1.001) / steps;
    let z = 0;
    
    for (let i = 0; i < steps; i++) {
        const ri = b0 * 1.001 + (i + 0.5) * dr;
        const bi = config.b_func(ri);
        if (ri > bi) {
            z += Math.sqrt(bi / (ri - bi)) * dr;
        }
    }
    return z;
}

function createWormholeGeometry() {
    const config = getShapeConfig();
    const b0 = params.b0;
    const rMin = b0, rMax = b0 * 6;
    const rSegments = 60, thetaSegments = 48;
    
    const geometry = new THREE.BufferGeometry();
    const vertices = [], normals = [], colors = [], indices = [];
    
    for (let side = -1; side <= 1; side += 2) {
        for (let i = 0; i <= rSegments; i++) {
            const t = i / rSegments;
            const r = rMin + t * (rMax - rMin);
            const z = side * embeddingZ(r, config);
            
            for (let j = 0; j <= thetaSegments; j++) {
                const theta = (j / thetaSegments) * Math.PI * 2;
                vertices.push(r * Math.cos(theta), z, r * Math.sin(theta));
                
                const dr = 0.01;
                const z1 = side * embeddingZ(r + dr, config);
                const dz_dr = (z1 - z) / dr;
                const nx = -dz_dr * Math.cos(theta), ny = 1, nz = -dz_dr * Math.sin(theta);
                const len = Math.sqrt(nx*nx + ny*ny + nz*nz);
                normals.push(nx/len, ny/len, nz/len);
                
                if (params.vizMode === 'curvature') {
                    const b = config.b_func(r);
                    const bp = config.b_deriv(r);
                    const K = Math.abs(bp * b / (2 * r * r * r));
                    const intensity = Math.min(K * 10, 1);
                    colors.push(0.3 + 0.7 * intensity, 0.5 - 0.3 * intensity, 1 - 0.5 * intensity);
                } else {
                    const dist = (r - b0) / (rMax - b0);
                    colors.push(0.3 + 0.4 * dist, 0.5 + 0.3 * dist, 0.9 - 0.3 * dist);
                }
            }
        }
    }
    
    const vertsPerRing = thetaSegments + 1;
    const vertsPerSide = (rSegments + 1) * vertsPerRing;
    
    for (let side = 0; side < 2; side++) {
        const offset = side * vertsPerSide;
        for (let i = 0; i < rSegments; i++) {
            for (let j = 0; j < thetaSegments; j++) {
                const a = offset + i * vertsPerRing + j;
                const b = offset + i * vertsPerRing + j + 1;
                const c = offset + (i + 1) * vertsPerRing + j;
                const d = offset + (i + 1) * vertsPerRing + j + 1;
                indices.push(a, c, b, b, c, d);
            }
        }
    }
    
    for (let j = 0; j < thetaSegments; j++) {
        indices.push(j, j + 1, vertsPerSide + j);
        indices.push(j + 1, vertsPerSide + j + 1, vertsPerSide + j);
    }
    
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
    geometry.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.setIndex(indices);
    
    return geometry;
}

function updateWormhole() {
    if (wormholeMesh) {
        scene.remove(wormholeMesh);
        wormholeMesh.geometry.dispose();
        wormholeMesh.material.dispose();
    }
    
    const geometry = createWormholeGeometry();
    const material = new THREE.MeshPhongMaterial({
        vertexColors: true,
        side: THREE.DoubleSide,
        shininess: 30,
        transparent: true,
        opacity: 0.9
    });
    wormholeMesh = new THREE.Mesh(geometry, material);
    scene.add(wormholeMesh);
    
    if (gridHelper) scene.remove(gridHelper);
    gridHelper = new THREE.GridHelper(20, 20, 0x333355, 0x222233);
    gridHelper.position.y = -5;
    scene.add(gridHelper);
    
    geodesicLines.forEach(line => scene.remove(line));
    geodesicLines = [];
    if (params.vizMode === 'geodesic') addGeodesics();
}

function addGeodesics() {
    const config = getShapeConfig();
    const b0 = params.b0;
    const colors = [0xff6666, 0x66ff66, 0x6666ff, 0xffff66];
    
    for (let g = 0; g < 4; g++) {
        const points = [];
        const startAngle = (g / 4) * Math.PI * 2;
        
        for (let i = 0; i <= 100; i++) {
            const t = i / 100;
            const r = b0 + t * (b0 * 5);
            const theta = startAngle + t * 0.5;
            points.push(new THREE.Vector3(r * Math.cos(theta), embeddingZ(r, config), r * Math.sin(theta)));
        }
        
        const line = new THREE.Line(
            new THREE.BufferGeometry().setFromPoints(points),
            new THREE.LineBasicMaterial({ color: colors[g] })
        );
        scene.add(line);
        geodesicLines.push(line);
    }
}

// ============================================================================
// RAYCASTING & CLICK HANDLING
// ============================================================================

const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();
let clickMarker = null;

function onCanvasClick(event) {
    const rect = canvas.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    
    raycaster.setFromCamera(mouse, camera);
    
    if (wormholeMesh) {
        const intersects = raycaster.intersectObject(wormholeMesh);
        if (intersects.length > 0) {
            const point = intersects[0].point;
            const r = Math.sqrt(point.x * point.x + point.z * point.z);
            const theta = Math.atan2(point.z, point.x);
            
            document.getElementById('click-coords').textContent = 
                `r = ${r.toFixed(2)}, θ = ${(theta * 180 / Math.PI).toFixed(0)}°`;
            
            displayDerivation(r, theta);
            addClickMarker(point);
        }
    }
}

function addClickMarker(position) {
    if (clickMarker) scene.remove(clickMarker);
    
    clickMarker = new THREE.Mesh(
        new THREE.SphereGeometry(0.1, 16, 16),
        new THREE.MeshBasicMaterial({ color: 0xff7700 })
    );
    clickMarker.position.copy(position);
    scene.add(clickMarker);
    
    setTimeout(() => {
        if (clickMarker) {
            scene.remove(clickMarker);
            clickMarker = null;
        }
    }, 3000);
}

// ============================================================================
// DERIVATION DISPLAY WITH KATEX
// ============================================================================

function renderLatex(latex, element) {
    try {
        katex.render(latex, element, { throwOnError: false, displayMode: false });
    } catch (e) {
        element.textContent = latex;
    }
}

function displayDerivation(r, theta) {
    const config = getShapeConfig();
    const b0 = params.b0;
    const b = config.b_func(r);
    const bp = config.b_deriv(r);
    
    const content = document.getElementById('derivation-content');
    content.innerHTML = '';
    
    // Section 1: Point & Configuration
    let html = `
        <div class="derivation-section">
            <h3>Point & Configuration</h3>
            <div class="params-grid">
                <div class="param-item">
                    <span class="label">Position</span>
                    <span class="value">r = ${r.toFixed(3)}</span>
                </div>
                <div class="param-item">
                    <span class="label">Angle</span>
                    <span class="value">θ = ${(theta * 180 / Math.PI).toFixed(1)}°</span>
                </div>
                <div class="param-item">
                    <span class="label">Shape</span>
                    <span class="value">${config.name}</span>
                </div>
                <div class="param-item">
                    <span class="label">Throat</span>
                    <span class="value">b₀ = ${b0}</span>
                </div>
            </div>
        </div>
    `;
    
    // Section 2: Metric tensor derivation
    const g_rr = r / (r - b);
    const g_thth = r * r;
    const g_phph = r * r * Math.sin(theta) * Math.sin(theta);
    
    html += `
        <div class="derivation-section">
            <h3>Metric Tensor</h3>
            <div class="derivation-line" id="metric-line-1"></div>
            <div class="derivation-line" id="metric-line-2"></div>
            <div class="derivation-line" id="metric-line-3"></div>
            <div class="derivation-line" id="metric-line-4"></div>
            <div class="derivation-line" id="metric-line-5"></div>
        </div>
    `;
    
    // Section 3: Christoffel derivation
    html += `
        <div class="derivation-section">
            <h3>Christoffel Symbol Γʳᵣᵣ</h3>
            <div class="derivation-line" id="chris-line-1"></div>
            <div class="derivation-line" id="chris-line-2"></div>
            <div class="derivation-line" id="chris-line-3"></div>
            <div class="derivation-line" id="chris-line-4"></div>
            <div class="derivation-line" id="chris-line-5"></div>
            <div class="derivation-line" id="chris-line-6"></div>
        </div>
    `;
    
    // Section 4: Numerical values
    const christoffels = computeChristoffelValues(r, theta, config);
    
    html += `
        <div class="derivation-section">
            <h3>All Christoffel Symbols</h3>
            <div class="tensor-grid">
                ${christoffels.map(c => `
                    <div class="tensor-item">
                        <span class="name">${c.name}</span>
                        <span class="val">${c.value.toFixed(4)}</span>
                    </div>
                `).join('')}
            </div>
        </div>
    `;
    
    // Section 5: Curvature
    const K = -bp * b / (2 * r * r * r);
    
    html += `
        <div class="derivation-section">
            <h3>Curvature</h3>
            <div class="derivation-line" id="curv-line-1"></div>
            <div class="derivation-line" id="curv-line-2"></div>
        </div>
    `;
    
    content.innerHTML = html;
    
    // Now render LaTeX into the placeholder divs
    setTimeout(() => {
        // Metric lines
        renderLatex(`ds^2 = -dt^2 + \\frac{r}{r - b(r)}dr^2 + r^2 d\\theta^2 + r^2\\sin^2\\theta\\, d\\phi^2`, 
            document.getElementById('metric-line-1'));
        renderLatex(`g_{tt} = -1, \\quad g_{rr} = \\frac{r}{r - b(r)} = \\frac{${r.toFixed(2)}}{${r.toFixed(2)} - ${b.toFixed(2)}} = ${g_rr.toFixed(4)}`, 
            document.getElementById('metric-line-2'));
        renderLatex(`g_{\\theta\\theta} = r^2 = ${g_thth.toFixed(4)}`, 
            document.getElementById('metric-line-3'));
        renderLatex(`g_{\\phi\\phi} = r^2\\sin^2\\theta = ${g_phph.toFixed(4)}`, 
            document.getElementById('metric-line-4'));
        renderLatex(`g^{rr} = \\frac{r - b(r)}{r} = ${(1/g_rr).toFixed(4)}`, 
            document.getElementById('metric-line-5'));
        
        // Christoffel derivation
        renderLatex(`\\Gamma^\\sigma_{\\mu\\nu} = \\frac{1}{2}g^{\\sigma\\rho}\\left(\\partial_\\mu g_{\\nu\\rho} + \\partial_\\nu g_{\\mu\\rho} - \\partial_\\rho g_{\\mu\\nu}\\right)`, 
            document.getElementById('chris-line-1'));
        renderLatex(`\\Gamma^r_{rr} = \\frac{1}{2}g^{rr}\\left(\\partial_r g_{rr} + \\partial_r g_{rr} - \\partial_r g_{rr}\\right) = \\frac{1}{2}g^{rr}\\partial_r g_{rr}`, 
            document.getElementById('chris-line-2'));
        renderLatex(`\\partial_r g_{rr} = \\partial_r\\left(\\frac{r}{r-b}\\right) = \\frac{(r-b) - r(1-b')}{(r-b)^2} = \\frac{rb' - b}{(r-b)^2}`, 
            document.getElementById('chris-line-3'));
        
        const dg_rr = (r * bp - b) / ((r - b) * (r - b));
        renderLatex(`\\partial_r g_{rr} = \\frac{${r.toFixed(2)} \\cdot ${bp.toFixed(4)} - ${b.toFixed(2)}}{(${r.toFixed(2)} - ${b.toFixed(2)})^2} = ${dg_rr.toFixed(4)}`, 
            document.getElementById('chris-line-4'));
        
        const gamma_r_rr = 0.5 * (1/g_rr) * dg_rr;
        renderLatex(`\\Gamma^r_{rr} = \\frac{1}{2} \\cdot ${(1/g_rr).toFixed(4)} \\cdot ${dg_rr.toFixed(4)} = ${gamma_r_rr.toFixed(4)}`, 
            document.getElementById('chris-line-5'));
        
        renderLatex(`\\boxed{\\Gamma^r_{rr} = \\frac{b - rb'}{2r(r-b)} = ${gamma_r_rr.toFixed(6)}}`, 
            document.getElementById('chris-line-6'));
        
        // Curvature
        renderLatex(`R = 0 \\quad \\text{(vacuum solution, Morris-Thorne)}`, 
            document.getElementById('curv-line-1'));
        renderLatex(`K_{\\text{Gaussian}} = -\\frac{b \\cdot b'}{2r^3} = ${K.toFixed(6)}`, 
            document.getElementById('curv-line-2'));
        
        // Animate highlighting
        document.querySelectorAll('.derivation-line').forEach((el, i) => {
            setTimeout(() => {
                el.classList.add('highlight');
                setTimeout(() => el.classList.remove('highlight'), 400);
            }, i * 150);
        });
        
    }, 50);
}

function computeChristoffelValues(r, theta, config) {
    const b = config.b_func(r);
    const bp = config.b_deriv(r);
    const sinT = Math.sin(theta);
    const cosT = Math.cos(theta);
    const rmb = r - b;
    
    const values = [];
    
    if (Math.abs(rmb) > 0.001) {
        values.push({ name: 'Γʳᵣᵣ', value: (b - r * bp) / (2 * r * rmb) });
    }
    values.push({ name: 'Γʳθθ', value: -(r - b) });
    values.push({ name: 'Γʳφφ', value: -(r - b) * sinT * sinT });
    values.push({ name: 'Γᶿᵣθ', value: 1 / r });
    values.push({ name: 'Γᶿφφ', value: -sinT * cosT });
    values.push({ name: 'Γᵠᵣφ', value: 1 / r });
    if (Math.abs(sinT) > 0.001) {
        values.push({ name: 'Γᵠθφ', value: cosT / sinT });
    }
    
    return values;
}

// ============================================================================
// CAMERA CONTROLS
// ============================================================================

let isDragging = false;
let previousMousePosition = { x: 0, y: 0 };
let cameraAngle = { theta: Math.PI / 4, phi: Math.PI / 4 };
let cameraDistance = 12;

function updateCameraPosition() {
    camera.position.x = cameraDistance * Math.sin(cameraAngle.theta) * Math.cos(cameraAngle.phi);
    camera.position.y = cameraDistance * Math.sin(cameraAngle.phi);
    camera.position.z = cameraDistance * Math.cos(cameraAngle.theta) * Math.cos(cameraAngle.phi);
    camera.lookAt(0, 0, 0);
}

canvas.addEventListener('mousedown', (e) => {
    if (e.button === 2 || e.shiftKey) {
        isDragging = true;
        previousMousePosition = { x: e.clientX, y: e.clientY };
    } else {
        onCanvasClick(e);
    }
});

canvas.addEventListener('mousemove', (e) => {
    if (isDragging) {
        const deltaX = e.clientX - previousMousePosition.x;
        const deltaY = e.clientY - previousMousePosition.y;
        cameraAngle.theta += deltaX * 0.01;
        cameraAngle.phi = Math.max(-Math.PI/2 + 0.1, Math.min(Math.PI/2 - 0.1, cameraAngle.phi + deltaY * 0.01));
        updateCameraPosition();
        previousMousePosition = { x: e.clientX, y: e.clientY };
    }
});

canvas.addEventListener('mouseup', () => isDragging = false);
canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    cameraDistance = Math.max(5, Math.min(30, cameraDistance + e.deltaY * 0.01));
    updateCameraPosition();
});
canvas.addEventListener('contextmenu', (e) => e.preventDefault());

// ============================================================================
// UI CONTROLS
// ============================================================================

document.getElementById('shape-select').addEventListener('change', (e) => {
    params.shapeType = e.target.value;
    updateWormhole();
});

document.getElementById('throat-radius').addEventListener('input', (e) => {
    params.b0 = parseFloat(e.target.value);
    document.getElementById('throat-value').textContent = params.b0.toFixed(1);
    updateWormhole();
});

document.getElementById('viz-mode').addEventListener('change', (e) => {
    params.vizMode = e.target.value;
    updateWormhole();
});

// ============================================================================
// ANIMATION & INIT
// ============================================================================

function animate() {
    requestAnimationFrame(animate);
    renderer.render(scene, camera);
}

window.addEventListener('resize', () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
});

updateCameraPosition();
updateWormhole();
animate();

console.log('Wormhole Tensor Calculator initialized');
