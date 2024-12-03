import React, { useEffect, useRef } from 'react';
import { Atoms, WEAS } from 'weas';

// Define the useNoDrag hook
function useNoDrag(ref, disabled = false) {
  useEffect(() => {
      const handleClick = (e) => {
          if (disabled) return;

          e.stopPropagation(); // Stop the event from bubbling up to prevent dragging
      };
      const el = ref.current;
      el?.addEventListener('pointerdown', handleClick);

      return () => {
          el?.removeEventListener('pointerdown', handleClick);
      };
  }, [ref, disabled]); // Re-run when ref or disabled changes
}

function AtomsItem({ data }) {
  const weasContainerRef = useRef(null);
  useNoDrag(weasContainerRef); // Apply the useNoDrag hook


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

    // console.log("atoms control data: ", data)
    data = data.data;
    let atomsData = {};
    if (data.node_type === 'data.core.structure.StructureData.') {
      atomsData = structureToAtomsData(data)
    } else if (data.node_type === 'data.workgraph.ase.atoms.Atoms.AtomsData.') {
      atomsData = aseAtomsToAtomsData(data)
    }
    // console.log("atoms control data: ", atomsData);
    const atoms = new Atoms(atomsData);

    if (weasContainerRef.current) {
      const defaultGuiConfig = {
        controls: {
           enabled: false,
           atomsControl: false,
           colorControl: false,
           cameraControls: false,
        },
        buttons: {
           enabled: true,
           fullscreen: true,
          //  undo: false,
          //  redo: false,
           download: false,
          //  measurement: false,
        }
        };

        function preventEventPropagation(element) {
          const stopPropagation = (e) => e.stopPropagation();
          ["click", "keydown", "keyup", "keypress"].forEach((eventType) => {
            element.addEventListener(eventType, stopPropagation, false);
          });
        }

        let domElement = document.createElement("div");
        domElement.style.cssText = "position: relative; width: 195px; height: 195px; border: 1px solid black;";
        weasContainerRef.current.appendChild(domElement);

      // Create an instance of AtomsViewer and pass the Atoms object to it
      preventEventPropagation(domElement);
      const editor = new WEAS({domElement: domElement, guiConfig: defaultGuiConfig});

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
      {/* <h1>Atoms Viewer</h1> */}
      <div ref={weasContainerRef} style={{ position: "relative", width: '200px', height: '200px' }}></div>
    </div>
  );
}

export default AtomsItem;
