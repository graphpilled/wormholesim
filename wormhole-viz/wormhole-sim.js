// wormhole-sim.js - Full wormhole simulator with geodesic physics
// Now uses OCaml TensorEngine compiled via js_of_ocaml for actual tensor computation

// ============================================================================
// INITIALIZATION
// ============================================================================

const canvas = document.getElementById('wormhole-canvas');
const container = document.getElementById('canvas-container');

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020208);

const camera = new THREE.PerspectiveCamera(60, container.clientWidth / container.clientHeight, 0.1, 1000);
camera.position.set(8, 5, 8);
camera.lookAt(0, 0, 0);

const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setSize(container.clientWidth, container.clientHeight);
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

// Lighting
const ambientLight = new THREE.AmbientLight(0x404060, 0.4);
scene.add(ambientLight);

const directionalLight = new THREE.DirectionalLight(0xffffff, 0.6);
directionalLight.position.set(5, 10, 5);
scene.add(directionalLight);

// Throat glow
const throatLight = new THREE.PointLight(0x4488ff, 1.5, 8);
throatLight.position.set(0, 0, 0);
scene.add(throatLight);

// ============================================================================
// SIMULATION STATE
// ============================================================================

let wormholeMesh = null;
let accretionDisk = null;
let photonSphere = null;
let gridHelpers = [];
let universeMarkers = [];

const particles = [];
const lightRays = [];
const trajectoryLines = [];

let simTime = 0;
let properTime = 0;
let coordTime = 0;
let throughThroatCount = 0;

let firstPersonMode = false;
let firstPersonParticle = null;
let savedCameraState = null;

const params = {
    shapeType: 'constant',
    b0: 1.0,
    vizMode: 'embedding',
    simSpeed: 1.0
};

// ============================================================================
// OCAML TENSOR ENGINE BRIDGE
// ============================================================================

// Check TensorEngine is loaded
if (typeof TensorEngine === 'undefined') {
    console.error('TensorEngine not loaded! Include tensor_engine.js before this script.');
}

// Wrapper to get Christoffel symbols from OCaml engine
function christoffel(r, theta, config) {
    const G = TensorEngine.computeChristoffel(r, theta, params.b0);
    // Map from camelCase (OCaml) to underscore format (existing code)
    return {
        r_rr: G.rRr,
        r_thth: G.rThth,
        r_phph: G.rPhph,
        th_rth: G.thRth,
        th_phph: G.thPhph,
        ph_rph: G.phRph,
        ph_thph: G.phThph
    };
}

// Wrapper for embedding function
function embeddingZ(r, config, side = 1) {
    return TensorEngine.embeddingZ(params.b0, r) * side;
}

// Wrapper for shape function  
function getShapeAt(r) {
    return TensorEngine.getShapeAt(params.b0, r);
}

// ============================================================================
// WORMHOLE PHYSICS (using TensorEngine)
// ============================================================================

function getShapeConfig() {
    const b0 = params.b0;
    switch (params.shapeType) {
        case 'ellis':
            return {
                name: `Ellis b(r) = ${b0}²/r`,
                b_func: (r) => b0 * b0 / r,
                b_deriv: (r) => -b0 * b0 / (r * r)
            };
        case 'schwarzschild':
            return {
                name: `Schwarzschild-like`,
                b_func: (r) => b0,
                b_deriv: (r) => 0,
                hasSingularity: true
            };
        default:
            return {
                name: `Morris-Thorne b(r) = ${b0}`,
                b_func: (r) => b0,
                b_deriv: (r) => 0
            };
    }
}

// Time dilation factor
function timeDilationFactor(r, config) {
    const b = getShapeAt(r);
    if (r <= b) return Infinity;
    return Math.sqrt(r / (r - b));
}

// ============================================================================
// GEODESIC INTEGRATION (using TensorEngine Christoffels)
// ============================================================================

// Geodesic equation RHS using OCaml-computed Christoffel symbols
function geodesicRHS(state, config, isNull = false) {
    const [t, r, theta, phi, ut, ur, uth, uph] = state;
    
    const rSafe = Math.max(r, params.b0 * 1.001);
    const G = christoffel(rSafe, theta, config);
    
    // dx^μ/dτ = u^μ
    const dt = ut;
    const dr = ur;
    const dth = uth;
    const dph = uph;
    
    // du^μ/dτ = -Γ^μ_αβ u^α u^β
    const dut = 0; // Φ = 0 means Γ^t_αβ = 0
    
    const dur = -G.r_rr * ur * ur 
                - G.r_thth * uth * uth 
                - G.r_phph * uph * uph;
    
    const duth = -2 * G.th_rth * ur * uth 
                 - G.th_phph * uph * uph;
    
    const duph = -2 * G.ph_rph * ur * uph 
                 - 2 * G.ph_thph * uth * uph;
    
    return [dt, dr, dth, dph, dut, dur, duth, duph];
}

// RK4 integration step
function rk4Step(state, dt, config, isNull = false) {
    const k1 = geodesicRHS(state, config, isNull);
    
    const state2 = state.map((v, i) => v + 0.5 * dt * k1[i]);
    const k2 = geodesicRHS(state2, config, isNull);
    
    const state3 = state.map((v, i) => v + 0.5 * dt * k2[i]);
    const k3 = geodesicRHS(state3, config, isNull);
    
    const state4 = state.map((v, i) => v + dt * k3[i]);
    const k4 = geodesicRHS(state4, config, isNull);
    
    return state.map((v, i) => 
        v + (dt / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
    );
}

// Normalize 4-velocity for timelike geodesic
function normalizeTimelike(state, config) {
    let [t, r, theta, phi, ut, ur, uth, uph] = state;
    
    const b = getShapeAt(r);
    const sinT = Math.sin(theta);
    
    const g_rr = r / (r - b);
    const g_thth = r * r;
    const g_phph = r * r * sinT * sinT;
    
    const spatial = g_rr * ur * ur + g_thth * uth * uth + g_phph * uph * uph;
    const utSq = 1 + spatial;
    
    if (utSq > 0) {
        ut = Math.sqrt(utSq);
    }
    
    return [t, r, theta, phi, ut, ur, uth, uph];
}

// Normalize for null geodesic
function normalizeNull(state, config) {
    let [t, r, theta, phi, kt, kr, kth, kph] = state;
    
    const b = getShapeAt(r);
    const sinT = Math.sin(theta);
    
    const g_rr = r / (r - b);
    const g_thth = r * r;
    const g_phph = r * r * sinT * sinT;
    
    const spatial = g_rr * kr * kr + g_thth * kth * kth + g_phph * kph * kph;
    kt = Math.sqrt(Math.max(0, spatial));
    
    return [t, r, theta, phi, kt, kr, kth, kph];
}

// ============================================================================
// PARTICLE & LIGHT RAY OBJECTS
// ============================================================================

class Particle {
    constructor(r, theta, phi, ur, uth, uph, color = 0xff4444) {
        this.config = getShapeConfig();
        
        this.state = [0, r, theta, phi, 0, ur, uth, uph];
        this.state = normalizeTimelike(this.state, this.config);
        
        this.side = 1;
        this.properTime = 0;
        this.alive = true;
        this.crossedThroat = false;
        
        this.trail = [];
        this.maxTrailLength = 200;
        
        this.mesh = new THREE.Mesh(
            new THREE.SphereGeometry(0.08, 16, 16),
            new THREE.MeshBasicMaterial({ color })
        );
        scene.add(this.mesh);
        
        this.trailGeometry = new THREE.BufferGeometry();
        this.trailLine = new THREE.Line(
            this.trailGeometry,
            new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.6 })
        );
        scene.add(this.trailLine);
        
        this.glow = new THREE.Mesh(
            new THREE.SphereGeometry(0.15, 8, 8),
            new THREE.MeshBasicMaterial({ color, transparent: true, opacity: 0.3 })
        );
        scene.add(this.glow);
        
        this.updateVisual();
    }
    
    step(dt) {
        if (!this.alive) return;
        
        const config = this.config;
        const b0 = params.b0;
        
        this.state = rk4Step(this.state, dt, config, false);
        this.properTime += dt;
        
        let [t, r, theta, phi, ut, ur, uth, uph] = this.state;
        
        if (r < b0 * 1.01 && !this.crossedThroat) {
            this.side *= -1;
            this.crossedThroat = true;
            throughThroatCount++;
            r = b0 * 1.01;
            ur = Math.abs(ur) * this.side * -1;
        } else if (r > b0 * 1.5) {
            this.crossedThroat = false;
        }
        
        if (theta < 0.01) { theta = 0.01; uth = Math.abs(uth); }
        if (theta > Math.PI - 0.01) { theta = Math.PI - 0.01; uth = -Math.abs(uth); }
        phi = ((phi % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
        
        this.state = [t, r, theta, phi, ut, ur, uth, uph];
        
        if (r > b0 * 8) {
            this.alive = false;
        }
        
        this.updateVisual();
    }
    
    updateVisual() {
        const [t, r, theta, phi] = this.state;
        const config = this.config;
        
        const z = embeddingZ(r, config, this.side);
        const x = r * Math.cos(phi);
        const y = z;
        const zPos = r * Math.sin(phi);
        
        this.mesh.position.set(x, y, zPos);
        this.glow.position.set(x, y, zPos);
        
        this.trail.push(new THREE.Vector3(x, y, zPos));
        if (this.trail.length > this.maxTrailLength) {
            this.trail.shift();
        }
        
        if (this.trail.length > 1) {
            this.trailGeometry.setFromPoints(this.trail);
        }
    }
    
    remove() {
        scene.remove(this.mesh);
        scene.remove(this.trailLine);
        scene.remove(this.glow);
        this.mesh.geometry.dispose();
        this.mesh.material.dispose();
        this.trailGeometry.dispose();
        this.trailLine.material.dispose();
        this.glow.geometry.dispose();
        this.glow.material.dispose();
    }
    
    getPosition() {
        const [t, r, theta, phi] = this.state;
        const z = embeddingZ(r, this.config, this.side);
        return new THREE.Vector3(r * Math.cos(phi), z, r * Math.sin(phi));
    }
    
    getVelocity() {
        const [t, r, theta, phi, ut, ur, uth, uph] = this.state;
        return Math.sqrt(ur*ur + r*r*uth*uth + r*r*Math.sin(theta)**2*uph*uph) / ut;
    }
}

class LightRay {
    constructor(r, theta, phi, kr, kth, kph, color = 0xffff44) {
        this.config = getShapeConfig();
        
        this.state = [0, r, theta, phi, 0, kr, kth, kph];
        this.state = normalizeNull(this.state, this.config);
        
        this.side = 1;
        this.alive = true;
        this.crossedThroat = false;
        
        this.trail = [];
        this.maxTrailLength = 300;
        
        this.trailGeometry = new THREE.BufferGeometry();
        this.trailLine = new THREE.Line(
            this.trailGeometry,
            new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.8, linewidth: 2 })
        );
        scene.add(this.trailLine);
        
        this.head = new THREE.Mesh(
            new THREE.SphereGeometry(0.04, 8, 8),
            new THREE.MeshBasicMaterial({ color })
        );
        scene.add(this.head);
        
        this.updateVisual();
    }
    
    step(dt) {
        if (!this.alive) return;
        
        const config = this.config;
        const b0 = params.b0;
        
        for (let i = 0; i < 3; i++) {
            this.state = rk4Step(this.state, dt, config, true);
            
            let [t, r, theta, phi, kt, kr, kth, kph] = this.state;
            
            if (r < b0 * 1.01 && !this.crossedThroat) {
                this.side *= -1;
                this.crossedThroat = true;
                throughThroatCount++;
                r = b0 * 1.01;
                kr = Math.abs(kr) * this.side * -1;
            } else if (r > b0 * 1.5) {
                this.crossedThroat = false;
            }
            
            if (theta < 0.01) { theta = 0.01; kth = Math.abs(kth); }
            if (theta > Math.PI - 0.01) { theta = Math.PI - 0.01; kth = -Math.abs(kth); }
            phi = ((phi % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
            
            this.state = [t, r, theta, phi, kt, kr, kth, kph];
            
            if (r > b0 * 10) {
                this.alive = false;
                break;
            }
        }
        
        this.updateVisual();
    }
    
    updateVisual() {
        const [t, r, theta, phi] = this.state;
        const config = this.config;
        
        const z = embeddingZ(r, config, this.side);
        const x = r * Math.cos(phi);
        const y = z;
        const zPos = r * Math.sin(phi);
        
        this.head.position.set(x, y, zPos);
        
        this.trail.push(new THREE.Vector3(x, y, zPos));
        if (this.trail.length > this.maxTrailLength) {
            this.trail.shift();
        }
        
        if (this.trail.length > 1) {
            this.trailGeometry.setFromPoints(this.trail);
        }
    }
    
    remove() {
        scene.remove(this.trailLine);
        scene.remove(this.head);
        this.trailGeometry.dispose();
        this.trailLine.material.dispose();
        this.head.geometry.dispose();
        this.head.material.dispose();
    }
}

// ============================================================================
// WORMHOLE GEOMETRY (using TensorEngine.embeddingZ)
// ============================================================================

function createWormholeGeometry() {
    const config = getShapeConfig();
    const b0 = params.b0;
    const rMin = b0, rMax = b0 * 6;
    const rSegments = 50, thetaSegments = 40;
    
    const geometry = new THREE.BufferGeometry();
    const vertices = [], normals = [], colors = [], indices = [];
    
    for (let side = -1; side <= 1; side += 2) {
        for (let i = 0; i <= rSegments; i++) {
            const t = i / rSegments;
            const r = rMin + t * (rMax - rMin);
            const z = embeddingZ(r, config, side);
            
            for (let j = 0; j <= thetaSegments; j++) {
                const theta = (j / thetaSegments) * Math.PI * 2;
                
                const x = r * Math.cos(theta);
                const y = z;
                const zPos = r * Math.sin(theta);
                
                vertices.push(x, y, zPos);
                
                // Normal
                const dr = 0.01;
                const z1 = embeddingZ(r + dr, config, side);
                const dz_dr = (z1 - z) / dr;
                
                const nx = 0;
                const ny = side;
                const nz_contrib = -dz_dr;
                const len = Math.sqrt(1 + dz_dr * dz_dr);
                normals.push(-dz_dr * Math.cos(theta) / len, side / len, -dz_dr * Math.sin(theta) / len);
                
                // Color based on curvature
                const b = getShapeAt(r);
                const curvature = Math.abs(b / (r * r));
                const intensity = Math.min(curvature * 3, 1);
                
                const cr = 0.1 + 0.2 * (1 - t) + 0.3 * intensity;
                const cg = 0.2 + 0.3 * (1 - t);
                const cb = 0.5 + 0.3 * t;
                
                colors.push(cr, cg, cb);
            }
        }
    }
    
    // Indices
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
    
    // Connect at throat
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
    // Clear deformations in tensor engine
    TensorEngine.clearDeformations();
    
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
        opacity: 0.85
    });
    
    wormholeMesh = new THREE.Mesh(geometry, material);
    scene.add(wormholeMesh);
    
    updateGrids();
}

function updateGrids() {
    gridHelpers.forEach(g => scene.remove(g));
    gridHelpers = [];
    
    const config = getShapeConfig();
    const zExtent = embeddingZ(params.b0 * 6, config, 1);
    
    // Upper universe
    const gridUp = new THREE.GridHelper(10, 10, 0x224444, 0x112222);
    gridUp.position.y = zExtent + 0.1;
    scene.add(gridUp);
    gridHelpers.push(gridUp);
    
    // Lower universe
    const gridDown = new THREE.GridHelper(10, 10, 0x442244, 0x221122);
    gridDown.position.y = -zExtent - 0.1;
    scene.add(gridDown);
    gridHelpers.push(gridDown);
}

// ============================================================================
// CAMERA
// ============================================================================

let isDragging = false;
let previousMousePosition = { x: 0, y: 0 };
let cameraDistance = 12;
let cameraAngle = { theta: Math.PI / 4, phi: Math.PI / 6 };

function updateCameraPosition() {
    camera.position.x = cameraDistance * Math.sin(cameraAngle.theta) * Math.cos(cameraAngle.phi);
    camera.position.y = cameraDistance * Math.sin(cameraAngle.phi);
    camera.position.z = cameraDistance * Math.cos(cameraAngle.theta) * Math.cos(cameraAngle.phi);
    camera.lookAt(0, 0, 0);
}

// ============================================================================
// ACTIONS
// ============================================================================

function dropParticle() {
    const r = params.b0 * 3 + Math.random() * params.b0 * 2;
    const theta = Math.PI / 2 + (Math.random() - 0.5) * 0.3;
    const phi = Math.random() * Math.PI * 2;
    const ur = -0.3 - Math.random() * 0.2;
    const uth = (Math.random() - 0.5) * 0.1;
    const uph = (Math.random() - 0.5) * 0.2;
    
    const colors = [0xff4444, 0x44ff44, 0x4444ff, 0xff44ff, 0xffff44, 0x44ffff];
    const color = colors[particles.length % colors.length];
    
    particles.push(new Particle(r, theta, phi, ur, uth, uph, color));
}

function shootLight() {
    const r = params.b0 * 5;
    const theta = Math.PI / 2;
    const phi = Math.random() * Math.PI * 2;
    const kr = -1;
    const kth = (Math.random() - 0.5) * 0.2;
    const kph = (Math.random() - 0.5) * 0.3;
    
    lightRays.push(new LightRay(r, theta, phi, kr, kth, kph));
}

function togglePhotonSphere() {
    if (photonSphere) {
        scene.remove(photonSphere);
        photonSphere.geometry.dispose();
        photonSphere.material.dispose();
        photonSphere = null;
    } else {
        const geometry = new THREE.SphereGeometry(params.b0 * 1.5, 32, 16);
        const material = new THREE.MeshBasicMaterial({
            color: 0xffff00,
            transparent: true,
            opacity: 0.1,
            wireframe: true
        });
        photonSphere = new THREE.Mesh(geometry, material);
        scene.add(photonSphere);
    }
}

function toggleAccretionDisk() {
    if (accretionDisk) {
        scene.remove(accretionDisk);
        accretionDisk.geometry.dispose();
        accretionDisk.material.dispose();
        accretionDisk = null;
    } else {
        const geometry = new THREE.RingGeometry(params.b0 * 1.5, params.b0 * 4, 64);
        const material = new THREE.ShaderMaterial({
            uniforms: {
                time: { value: 0 },
                innerRadius: { value: params.b0 * 1.5 },
                outerRadius: { value: params.b0 * 4 }
            },
            vertexShader: `
                varying vec2 vUv;
                void main() {
                    vUv = uv;
                    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
                }
            `,
            fragmentShader: `
                uniform float time;
                uniform float innerRadius;
                uniform float outerRadius;
                varying vec2 vUv;
                void main() {
                    float r = mix(innerRadius, outerRadius, vUv.x);
                    float temp = 1.0 / (r * r);
                    float brightness = temp * (0.5 + 0.5 * sin(vUv.y * 20.0 + time * 3.0));
                    vec3 color = vec3(1.0, 0.6, 0.2) * brightness;
                    gl_FragColor = vec4(color, brightness * 0.7);
                }
            `,
            transparent: true,
            side: THREE.DoubleSide
        });
        accretionDisk = new THREE.Mesh(geometry, material);
        accretionDisk.rotation.x = Math.PI / 2;
        scene.add(accretionDisk);
    }
}

function toggleFirstPerson() {
    if (firstPersonMode) {
        exitFirstPerson();
    } else if (particles.length > 0) {
        firstPersonMode = true;
        firstPersonParticle = particles[particles.length - 1];
        savedCameraState = {
            distance: cameraDistance,
            angle: { ...cameraAngle }
        };
    }
}

function exitFirstPerson() {
    firstPersonMode = false;
    firstPersonParticle = null;
    if (savedCameraState) {
        cameraDistance = savedCameraState.distance;
        cameraAngle = savedCameraState.angle;
        updateCameraPosition();
    }
}

function updateFirstPersonCamera() {
    if (!firstPersonParticle || !firstPersonParticle.alive) {
        exitFirstPerson();
        return;
    }
    
    const pos = firstPersonParticle.getPosition();
    const [t, r, theta, phi, ut, ur, uth, uph] = firstPersonParticle.state;
    
    // Look in direction of motion
    const dir = new THREE.Vector3(
        ur * Math.cos(phi) - r * uph * Math.sin(phi),
        0,
        ur * Math.sin(phi) + r * uph * Math.cos(phi)
    ).normalize();
    
    camera.position.copy(pos);
    camera.lookAt(pos.clone().add(dir));
}

function clearAll() {
    particles.forEach(p => p.remove());
    particles.length = 0;
    lightRays.forEach(l => l.remove());
    lightRays.length = 0;
    throughThroatCount = 0;
    properTime = 0;
    coordTime = 0;
}

function onCanvasClick(e) {
    dropParticle();
}

// ============================================================================
// DERIVATION PANEL (using TensorEngine.getDerivationSteps)
// ============================================================================

function renderLatex(latex, element) {
    if (typeof katex !== 'undefined' && element) {
        katex.render(latex, element, { throwOnError: false, displayMode: true });
    }
}

function showDerivation() {
    const config = getShapeConfig();
    const r = params.b0 * 2;
    const theta = Math.PI / 2;
    
    // Get derivation data from OCaml engine
    const D = TensorEngine.getDerivationSteps(r, theta, params.b0);
    
    // Get Christoffels
    const G = TensorEngine.computeChristoffel(r, theta, params.b0);
    
    const dilation = D.timeDilation;
    
    const panel = document.getElementById('derivation-panel');
    if (!panel) return;
    
    panel.innerHTML = `
        <h2>Tensor Derivation at r = ${r.toFixed(2)}</h2>
        <p class="engine-note">Computed by OCaml TensorEngine (js_of_ocaml)</p>
        
        <div class="derivation-section">
            <h3>Morris-Thorne Metric</h3>
            <div class="derivation-line" id="metric-eq"></div>
            <p>Shape function: b(r) = ${D.b.toFixed(4)}, b'(r) = ${D.bp.toFixed(6)}</p>
        </div>
        
        <div class="derivation-section">
            <h3>Christoffel Symbols Γᵅᵦᵧ</h3>
            <div class="derivation-line" id="chris-eq"></div>
            <div class="tensor-grid">
                <div class="tensor-item">
                    <span class="name">Γʳᵣᵣ</span>
                    <span class="val">${G.rRr.toFixed(4)}</span>
                </div>
                <div class="tensor-item">
                    <span class="name">Γʳθθ</span>
                    <span class="val">${G.rThth.toFixed(4)}</span>
                </div>
                <div class="tensor-item">
                    <span class="name">Γʳφφ</span>
                    <span class="val">${G.rPhph.toFixed(4)}</span>
                </div>
                <div class="tensor-item">
                    <span class="name">Γᶿᵣθ</span>
                    <span class="val">${G.thRth.toFixed(4)}</span>
                </div>
                <div class="tensor-item">
                    <span class="name">Γᶿφφ</span>
                    <span class="val">${G.thPhph.toFixed(4)}</span>
                </div>
                <div class="tensor-item">
                    <span class="name">Γᵠᵣφ</span>
                    <span class="val">${G.phRph.toFixed(4)}</span>
                </div>
            </div>
        </div>
        
        <div class="derivation-section">
            <h3>Curvature Invariants</h3>
            <div class="tensor-grid">
                <div class="tensor-item">
                    <span class="name">Gaussian K</span>
                    <span class="val">${D.gaussianK.toFixed(6)}</span>
                </div>
                <div class="tensor-item">
                    <span class="name">Ricci R</span>
                    <span class="val">${D.ricciR.toFixed(6)}</span>
                </div>
            </div>
        </div>
        
        <div class="derivation-section">
            <h3>Geodesic Equation</h3>
            <div class="derivation-line" id="geo-eq"></div>
            <p style="font-size: 10px; color: #555; margin-top: 8px;">
                Particles follow geodesics — computed via RK4 integration using the Christoffel symbols above.
            </p>
        </div>
        
        <div class="derivation-section">
            <h3>Time Dilation</h3>
            <div class="time-bar">
                <div class="time-bar-label">Proper time flows ${dilation.toFixed(2)}x slower here</div>
                <div class="time-bar-container">
                    <div class="time-bar-fill" style="width: ${Math.min(dilation/3 * 100, 100)}%"></div>
                </div>
            </div>
        </div>
    `;
    
    setTimeout(() => {
        renderLatex(`ds^2 = -dt^2 + \\frac{r}{r-b(r)}dr^2 + r^2 d\\Omega^2`, 
            document.getElementById('metric-eq'));
        renderLatex(`\\Gamma^\\sigma_{\\mu\\nu} = \\frac{1}{2}g^{\\sigma\\rho}(\\partial_\\mu g_{\\nu\\rho} + \\partial_\\nu g_{\\mu\\rho} - \\partial_\\rho g_{\\mu\\nu})`, 
            document.getElementById('chris-eq'));
        renderLatex(`\\frac{d^2 x^\\mu}{d\\tau^2} + \\Gamma^\\mu_{\\alpha\\beta}\\frac{dx^\\alpha}{d\\tau}\\frac{dx^\\beta}{d\\tau} = 0`, 
            document.getElementById('geo-eq'));
    }, 50);
}

// ============================================================================
// SIMULATION LOOP
// ============================================================================

function updateStatus() {
    const particleCount = document.getElementById('particle-count');
    const lightCount = document.getElementById('light-count');
    const properTimeEl = document.getElementById('proper-time');
    const coordTimeEl = document.getElementById('coord-time');
    const throatCountEl = document.getElementById('throat-count');
    
    if (particleCount) particleCount.textContent = particles.filter(p => p.alive).length;
    if (lightCount) lightCount.textContent = lightRays.filter(l => l.alive).length;
    if (properTimeEl) properTimeEl.textContent = properTime.toFixed(2);
    if (coordTimeEl) coordTimeEl.textContent = coordTime.toFixed(2);
    if (throatCountEl) throatCountEl.textContent = throughThroatCount;
    
    const alive = particles.filter(p => p.alive);
    if (alive.length > 0) {
        const config = getShapeConfig();
        const avgDilation = alive.reduce((sum, p) => {
            return sum + timeDilationFactor(p.state[1], config);
        }, 0) / alive.length;
        const dilElem = document.getElementById('dilation-factor');
        if (dilElem) {
            dilElem.textContent = avgDilation.toFixed(2) + 'x';
            if (avgDilation > 2) dilElem.className = 'value warning';
            else if (avgDilation > 5) dilElem.className = 'value danger';
            else dilElem.className = 'value';
        }
    }
}

let lastTime = performance.now();

function simulate() {
    const now = performance.now();
    const dt = Math.min((now - lastTime) / 1000, 0.05) * params.simSpeed;
    lastTime = now;
    
    simTime += dt;
    properTime += dt;
    coordTime += dt;
    
    particles.forEach(p => p.step(dt * 0.5));
    lightRays.forEach(l => l.step(dt * 0.3));
    
    for (let i = particles.length - 1; i >= 0; i--) {
        if (!particles[i].alive) {
            particles[i].remove();
            particles.splice(i, 1);
        }
    }
    
    for (let i = lightRays.length - 1; i >= 0; i--) {
        if (!lightRays[i].alive) {
            lightRays[i].remove();
            lightRays.splice(i, 1);
        }
    }
    
    if (accretionDisk) {
        accretionDisk.material.uniforms.time.value = simTime;
        accretionDisk.rotation.z = simTime * 0.2;
    }
    
    if (firstPersonMode) {
        updateFirstPersonCamera();
    }
    
    updateStatus();
}

function animate() {
    requestAnimationFrame(animate);
    simulate();
    renderer.render(scene, camera);
}

// ============================================================================
// EVENT LISTENERS
// ============================================================================

canvas.addEventListener('mousedown', (e) => {
    if (e.button === 2 || e.shiftKey) {
        isDragging = true;
        previousMousePosition = { x: e.clientX, y: e.clientY };
    } else {
        onCanvasClick(e);
    }
});

canvas.addEventListener('mousemove', (e) => {
    if (isDragging && !firstPersonMode) {
        const deltaX = e.clientX - previousMousePosition.x;
        const deltaY = e.clientY - previousMousePosition.y;
        cameraAngle.theta += deltaX * 0.01;
        cameraAngle.phi = Math.max(-Math.PI/2 + 0.1, Math.min(Math.PI/2 - 0.1, cameraAngle.phi + deltaY * 0.01));
        updateCameraPosition();
        previousMousePosition = { x: e.clientX, y: e.clientY };
    }
});

canvas.addEventListener('mouseup', () => isDragging = false);
canvas.addEventListener('mouseleave', () => isDragging = false);

canvas.addEventListener('wheel', (e) => {
    if (firstPersonMode) return;
    e.preventDefault();
    cameraDistance = Math.max(5, Math.min(30, cameraDistance + e.deltaY * 0.02));
    updateCameraPosition();
});

canvas.addEventListener('contextmenu', (e) => e.preventDefault());

// UI Controls
const shapeSelect = document.getElementById('shape-select');
if (shapeSelect) {
    shapeSelect.addEventListener('change', (e) => {
        params.shapeType = e.target.value;
        updateWormhole();
        particles.forEach(p => p.config = getShapeConfig());
        lightRays.forEach(l => l.config = getShapeConfig());
    });
}

const throatRadius = document.getElementById('throat-radius');
if (throatRadius) {
    throatRadius.addEventListener('input', (e) => {
        params.b0 = parseFloat(e.target.value);
        const throatValue = document.getElementById('throat-value');
        if (throatValue) throatValue.textContent = params.b0.toFixed(1);
        updateWormhole();
    });
}

const simSpeed = document.getElementById('sim-speed');
if (simSpeed) {
    simSpeed.addEventListener('input', (e) => {
        params.simSpeed = parseFloat(e.target.value);
        const speedValue = document.getElementById('speed-value');
        if (speedValue) speedValue.textContent = params.simSpeed.toFixed(1) + 'x';
    });
}

const vizMode = document.getElementById('viz-mode');
if (vizMode) {
    vizMode.addEventListener('change', (e) => {
        params.vizMode = e.target.value;
        updateWormhole();
    });
}

// Button listeners
const btnDropParticle = document.getElementById('btn-drop-particle');
if (btnDropParticle) btnDropParticle.addEventListener('click', dropParticle);

const btnShootLight = document.getElementById('btn-shoot-light');
if (btnShootLight) btnShootLight.addEventListener('click', shootLight);

const btnLightSphere = document.getElementById('btn-light-sphere');
if (btnLightSphere) btnLightSphere.addEventListener('click', togglePhotonSphere);

const btnFirstPerson = document.getElementById('btn-first-person');
if (btnFirstPerson) btnFirstPerson.addEventListener('click', toggleFirstPerson);

const btnAccretion = document.getElementById('btn-accretion');
if (btnAccretion) btnAccretion.addEventListener('click', toggleAccretionDisk);

const btnClear = document.getElementById('btn-clear');
if (btnClear) btnClear.addEventListener('click', clearAll);

// Keyboard
document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && firstPersonMode) {
        exitFirstPerson();
    }
    if (e.key === ' ') {
        e.preventDefault();
        dropParticle();
    }
    if (e.key === 'l') {
        shootLight();
    }
    if (e.key === 'f') {
        toggleFirstPerson();
    }
});

window.addEventListener('resize', () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
});

// ============================================================================
// INIT
// ============================================================================

updateCameraPosition();
updateWormhole();
animate();

console.log('🌀 Wormhole Simulator initialized');
console.log('📐 Using OCaml TensorEngine for geodesic computation');
console.log('Controls: Click = drop particle, Space = drop particle, L = light ray, F = first person, Esc = exit');
