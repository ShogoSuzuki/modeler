// standard global variables
var container, scene, camera, renderer, controls, stats;
var mesh, bbox_mesh, bbox_geo, bbox_material;
var bbox, x_max, x_min, x_range;
var y_max, y_min, y_range;
var z_max, z_min, z_range;


window.onload = main;


// The main function: init everything and animate everything
function main() {
    init();
    animate();
}


// The grid resolution: res * res * res
var grid_resolution = 64;


// Init a regular grid, sample the model,
// Generate a mesh using the Marching Cubes algorithm 
// for the sampled model.
// And set the ThreeJS renderer to render the mesh.
function init() {
    // Domain bounding box
    bbox = model.boundingBox();
    x_min = bbox["x_min"];
    x_max = bbox["x_max"];
    y_min = bbox["y_min"];
    y_max = bbox["y_max"];
    z_min = bbox["z_min"];
    z_max = bbox["z_max"];

    x_range = x_max - x_min;
    y_range = y_max - y_min;
    z_range = z_max - z_min;


    // SCENE
    scene = new THREE.Scene();

    // CAMERA
    var SCREEN_WIDTH = window.innerWidth, SCREEN_HEIGHT = window.innerHeight;
    var VIEW_ANGLE = 45, ASPECT = SCREEN_WIDTH / SCREEN_HEIGHT, NEAR = 0.1, FAR = 20000;
    camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
    scene.add(camera);
    camera.position.set(x_max + 0.5*x_range, 
			y_max + 0.5*y_range, 
			z_max + 0.5*z_range);
    camera.lookAt(scene.position);	


    // RENDERER
    if (Detector.webgl) {
	renderer = new THREE.WebGLRenderer({antialias:true});
    } else {
	renderer = new THREE.CanvasRenderer(); 
    }
    renderer.setSize(SCREEN_WIDTH, SCREEN_HEIGHT);
    renderer.setClearColor(0xffffff);
    container = document.getElementById('ThreeJS');
    container.appendChild(renderer.domElement);


    // EVENTS
    THREEx.WindowResize(renderer, camera);

    
    // CONTROLS
    controls = new THREE.OrbitControls( camera, renderer.domElement );


    // STATS
    stats = new Stats();
    stats.domElement.style.position = 'absolute';
    stats.domElement.style.bottom = '0px';
    stats.domElement.style.zIndex = 100;
    container.appendChild(stats.domElement);


    // lights
    // TODO: light position should depend on the scene geometry
    // for example we can place them at the corners of the bounding box
    var light1 = new THREE.PointLight(0xffffff);
    light1.position.set(0,y_max,0);
    scene.add(light1);

    var light2 = new THREE.PointLight(0xffffff);
    light2.position.set(x_max,0,0);
    scene.add(light2);

    var light3 = new THREE.PointLight(0xffffff);
    light3.position.set(0,0,z_max);
    scene.add(light3);

    var light4 = new THREE.PointLight(0xffffff);
    light4.position.set(0,y_min,0);
    scene.add(light4);

    var light5 = new THREE.PointLight(0xffffff);
    light5.position.set(x_min,0,0);
    scene.add(light5);

    var light6 = new THREE.PointLight(0xffffff);
    light6.position.set(0,0,x_min);
    scene.add(light6);
        
    scene.add( new THREE.AxisHelper(100) );


    // Add the polygonized implicit to the scene
    mesh = MC.polygonize(model, grid_resolution);
    scene.add(mesh);


    // add a grid to help position each object
    var grid = new THREE.GridHelper(500, 25);
    scene.add(grid);


    // Add the bounding box to the object
    bbox_geo = new THREE.BoxGeometry(x_range, y_range, z_range);
    bbox_material = new THREE.MeshBasicMaterial({color: 0x0000ff, wireframe: true});
    bbox_mesh = new THREE.Mesh(bbox_geo, bbox_material);
    scene.add(bbox_mesh);
}

function animate() 
{
    requestAnimationFrame( animate );
    render();		
    update();
}

function update()
{
    controls.update();
    stats.update();
}

function render() 
{
    renderer.render( scene, camera );
}
