var MC = {};


// Input: 
// * model: solid to be rendered (implicitly defined)
// * grid_resolution: number of subdivisions in each direction
//
// Output:
// * mesh object (threejs mesh)
//
MC.polygonize = function(model, grid_resolution) {

    var points = [];
    var values = [];


    // Domain bounding box
    var bbox = model.boundingBox();
    var x_min = bbox["x_min"];
    var x_max = bbox["x_max"];
    var y_min = bbox["y_min"];
    var y_max = bbox["y_max"];
    var z_min = bbox["z_min"];
    var z_max = bbox["z_max"];

    var x_range = x_max - x_min;
    var y_range = y_max - y_min;
    var z_range = z_max - z_min;


    // Grid resolution
    var size = grid_resolution;
    
    // Generate a list of 3D points and values at those points
    for (var k = 0; k < size; k++)
	for (var j = 0; j < size; j++)
	    for (var i = 0; i < size; i++)
    {
	// actual values
	var x = x_min + x_range * i / (size - 1);
	var y = y_min + y_range * j / (size - 1);
	var z = z_min + z_range * k / (size - 1);
	points.push( new THREE.Vector3(x,y,z) );
	var value = model.eval(x, y, z);
	values.push( value );
    }
    
    // Marching Cubes Algorithm
    
    var size2 = size * size;

    // Vertices may occur along edges of cube, when the values at the edge's endpoints
    //   straddle the isolevel value.
    // Actual position along edge weighted according to function values.
    var vlist = new Array(12);
    
    var geometry = new THREE.Geometry();
    var vertexIndex = 0;
    
    for (var z = 0; z < size - 1; z++)
	for (var y = 0; y < size - 1; y++)
	    for (var x = 0; x < size - 1; x++)
    {
	// index of base point, and also adjacent points on cube
	var p    = x + size * y + size2 * z,
	px   = p   + 1,
	py   = p   + size,
	pxy  = py  + 1,
	pz   = p   + size2,
	pxz  = px  + size2,
	pyz  = py  + size2,
	pxyz = pxy + size2;
	
	// store scalar values corresponding to vertices
	var value0 = values[ p    ],
	value1 = values[ px   ],
	value2 = values[ py   ],
	value3 = values[ pxy  ],
	value4 = values[ pz   ],
	value5 = values[ pxz  ],
	value6 = values[ pyz  ],
	value7 = values[ pxyz ];
	
	// place a "1" in bit positions corresponding to vertices whose
	//   isovalue is less than given constant.
	
	var isolevel = 0;
	
	var cubeindex = 0;
	if ( value0 < isolevel ) cubeindex |= 1;
	if ( value1 < isolevel ) cubeindex |= 2;
	if ( value2 < isolevel ) cubeindex |= 8;
	if ( value3 < isolevel ) cubeindex |= 4;
	if ( value4 < isolevel ) cubeindex |= 16;
	if ( value5 < isolevel ) cubeindex |= 32;
	if ( value6 < isolevel ) cubeindex |= 128;
	if ( value7 < isolevel ) cubeindex |= 64;
	
	// bits = 12 bit number, indicates which edges are crossed by the isosurface
	var bits = MC.edgeTable[ cubeindex ];
	
	// if none are crossed, proceed to next iteration
	if ( bits === 0 ) continue;
	
	// check which edges are crossed, and estimate the point location
	//    using a weighted average of scalar values at edge endpoints.
	// store the vertex in an array for use later.
	var mu = 0.5; 
	
	// bottom of the cube
	if ( bits & 1 )
	{		
	    mu = ( isolevel - value0 ) / ( value1 - value0 );
	    vlist[0] = points[p].clone().lerp( points[px], mu );
	}
	if ( bits & 2 )
	{
	    mu = ( isolevel - value1 ) / ( value3 - value1 );
	    vlist[1] = points[px].clone().lerp( points[pxy], mu );
	}
	if ( bits & 4 )
	{
	    mu = ( isolevel - value2 ) / ( value3 - value2 );
	    vlist[2] = points[py].clone().lerp( points[pxy], mu );
	}
	if ( bits & 8 )
	{
	    mu = ( isolevel - value0 ) / ( value2 - value0 );
	    vlist[3] = points[p].clone().lerp( points[py], mu );
	}
	// top of the cube
	if ( bits & 16 )
	{
	    mu = ( isolevel - value4 ) / ( value5 - value4 );
	    vlist[4] = points[pz].clone().lerp( points[pxz], mu );
	}
	if ( bits & 32 )
	{
	    mu = ( isolevel - value5 ) / ( value7 - value5 );
	    vlist[5] = points[pxz].clone().lerp( points[pxyz], mu );
	}
	if ( bits & 64 )
	{
	    mu = ( isolevel - value6 ) / ( value7 - value6 );
	    vlist[6] = points[pyz].clone().lerp( points[pxyz], mu );
	}
	if ( bits & 128 )
	{
	    mu = ( isolevel - value4 ) / ( value6 - value4 );
	    vlist[7] = points[pz].clone().lerp( points[pyz], mu );
	}
	// vertical lines of the cube
	if ( bits & 256 )
	{
	    mu = ( isolevel - value0 ) / ( value4 - value0 );
	    vlist[8] = points[p].clone().lerp( points[pz], mu );
	}
	if ( bits & 512 )
	{
	    mu = ( isolevel - value1 ) / ( value5 - value1 );
	    vlist[9] = points[px].clone().lerp( points[pxz], mu );
	}
	if ( bits & 1024 )
	{
	    mu = ( isolevel - value3 ) / ( value7 - value3 );
	    vlist[10] = points[pxy].clone().lerp( points[pxyz], mu );
	}
	if ( bits & 2048 )
	{
	    mu = ( isolevel - value2 ) / ( value6 - value2 );
	    vlist[11] = points[py].clone().lerp( points[pyz], mu );
	}
	
	// construct triangles -- get correct vertices from triTable.
	var i = 0;
	cubeindex <<= 4;  // multiply by 16... 
	// "Re-purpose cubeindex into an offset into triTable." 
	//  since each row really isn't a row.
	
	// the while loop should run at most 5 times,
	//   since the 16th entry in each row is a -1.
	while ( MC.triTable[ cubeindex + i ] != -1 ) 
	{
	    var index1 = MC.triTable[cubeindex + i];
	    var index2 = MC.triTable[cubeindex + i + 1];
	    var index3 = MC.triTable[cubeindex + i + 2];
	    
	    geometry.vertices.push( vlist[index1].clone() );
	    geometry.vertices.push( vlist[index2].clone() );
	    geometry.vertices.push( vlist[index3].clone() );
	    var face = 
		new THREE.Face3(vertexIndex, vertexIndex+1, vertexIndex+2);
	    geometry.faces.push( face );

	    vertexIndex += 3;
	    i += 3;
	}
    }
    
    geometry.computeFaceNormals();
    geometry.computeVertexNormals();
    
    var colorMaterial =  new THREE.MeshLambertMaterial( 
	{color: 0x7f7f7f, 
	 side: THREE.FrontSide,
	 shading: THREE.SmoothShading,
	 wireframe: false} );
    var mesh = new THREE.Mesh( geometry, colorMaterial );
    
    return mesh;
}
