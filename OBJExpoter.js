THREE.OBJExporter = function () {};

THREE.OBJExporter.prototype = {
    constructor: THREE.OBJExporter,

    parse: function ( object ) {

	var output = '';

	var vertex = new THREE.Vector3();

	var i, l;

        var mesh = object.mesh;

	var parseMesh = function ( mesh ) {
            var geometry = mesh.geometry;
            var settings = mesh.settings;

	    if ( geometry instanceof THREE.Geometry ) {

		// shortcuts
		var vertices = geometry.vertices;
		var indices = geometry.faces;

		// name of the mesh object
		output += 'o ' + mesh.name + '\n';

		// vertices

		if( vertices !== undefined ) {
		    for ( i = 0, l = vertices.length; i < l; i++ ) {
			vertex = vertices[i];

			// transform the vertex to export format
			output += 'v ' + vertex.x + ' ' + vertex.y + ' ' + vertex.z + '\n';
		    }
		}

		// faces
		if( indices !== null ) {
		    for ( i = 0, l = indices.length; i < l; i++ ) {
			output += 'f ' + (indices[i].a+1) + ' ' + (indices[i].b+1) + ' ' + (indices[i].c+1) + ' ' + "\n";
		    }
		} 
	    } else {
		console.warn( 'THREE.OBJExporter.parseMesh(): geometry type unsupported', geometry );
	    }
	};

	object.traverse( function ( child ) {
	    if ( child instanceof THREE.Mesh ) {
		parseMesh( child );
	    }
	} );

	return output;

    }

};

// Use FileSaver.js 'saveAs' function to save the string
function saveOBJ( scene, name ){
    var exporter = new THREE.OBJExporter();
    var objString = exporter.parse( scene );
    
    var blob = new Blob([objString], {type: 'text/plain'});
    

    function saveTextAsFile(blob, name) {
        var textToSaveAsURL = window.URL.createObjectURL(blob);
	
        var downloadLink = document.createElement("a");
        downloadLink.download = name;
        downloadLink.innerHTML = "Download File";
        downloadLink.href = textToSaveAsURL;
        downloadLink.onclick = destroyClickedElement;
        downloadLink.style.display = "none";
        document.body.appendChild(downloadLink);
	
        downloadLink.click();
    }

    function destroyClickedElement(event) {
        document.body.removeChild(event.target);
    }
    saveTextAsFile(blob, name + '.obj');
}