import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import AtomsItem from './AtomsItem.js'; // Adjust the path as necessary

import './DataNodeItem.css';
import '../App.css';

function DataNodeItem() {
  const { pk } = useParams();
  const [NodeData, setNodeData] = useState({ node_type: "" });

  useEffect(() => {
    fetch(`http://localhost:8000/api/datanode/${pk}`)
      .then(response => response.json())
      .then(data => {
        setNodeData(data);
      })
      .catch(error => console.error('Error fetching data:', error));
  }, [pk]); // Only re-run when `pk` changes

  // Safely convert any value to a string
  const stringifyValue = (value: unknown): string => {
    if (value === null) return 'null';
    if (typeof value === 'object') return JSON.stringify(value);
    return String(value);
  };

  return (
    <div className="table-container">
      <h2>DataNode</h2>
      <table>
        <thead>
          <tr>
            <th>Key</th>
            <th>Value</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(NodeData).map(([key, value]) => (
            <tr key={key}>
              <td>{key}</td>
              <td>{stringifyValue(value)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      {NodeData.node_type === 'data.core.structure.StructureData.' && <AtomsItem data={NodeData} />}
      {NodeData.node_type === 'data.workgraph.ase.atoms.Atoms.AtomsData.' && <AtomsItem data={NodeData} />}
    </div>
  );
}

export default DataNodeItem;
