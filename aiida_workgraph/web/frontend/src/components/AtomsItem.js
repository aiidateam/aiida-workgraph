import React, { useEffect, useRef } from 'react';
import { Atoms, WEAS } from 'weas';

function AtomsItem({ data }) {
  const weasContainerRef = useRef(null);


  // Convert AiiDA structure data to the format expected by Atoms
  function structureToAtomsData(inputData) {
    const data = {
      cell: inputData.cell,
      pbc: [inputData.pbc1, inputData.pbc2, inputData.pbc3],
      species: {},
      symbols: [],
      positions: []
    };

    // Process kinds to fill species information
    inputData.kinds.forEach((kind, index) => {
      data.species[kind.name] = kind.symbols[0]; // Assuming atomic number is the index + 1 for simplicity
    });

    // Process sites to fill positions and symbols
    inputData.sites.forEach(site => {
      data.symbols.push(site.kind_name); // Using index + 1 as a stand-in for atomic number
      data.positions.push(site.position);
    });

    return data;
}


  // Convert AiiDA structure data to the format expected by Atoms
  function aseAtomsToAtomsData(inputData) {
    const data = {
      cell: inputData.cell,
      pbc: inputData.pbc,
      symbols: inputData.symbols,
      positions: inputData.positions
    };

    return data;
}

  useEffect(() => {

    console.log("data: ", data)
    let atomsData = {};
    if (data.node_type === 'data.core.structure.StructureData.') {
      atomsData = structureToAtomsData(data)
    } else if (data.node_type === 'data.workgraph.ase.atoms.Atoms.AtomsData.') {
      atomsData = aseAtomsToAtomsData(data)
    }
    const atoms = new Atoms(atomsData);

    if (weasContainerRef.current) {
      // Create an instance of AtomsViewer and pass the Atoms object to it
      const editor = new WEAS({domElement: weasContainerRef.current});
      editor.avr.atoms = atoms;
      editor.render();

      // Cleanup function to be called when the component unmounts
      return () => {
        // viewer.destroy();
      };
    }
  }, [data]); // Include data in the dependency array

  return (
    <div>
      <h1>Atoms Viewer</h1>
      <div ref={weasContainerRef} style={{ position: "relative", width: '600px', height: '600px' }}></div>
    </div>
  );
}

export default AtomsItem;
