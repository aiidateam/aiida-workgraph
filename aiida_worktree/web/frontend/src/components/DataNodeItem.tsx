import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import './DataNodeItem.css';
import '../App.css';

function DataNodeItem() {
  const { pk } = useParams();
  const [NodeData, setNodeData] = useState({ summary: {}, nodes: {}, links: [], logs: [], pk: [] });

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
    </div>
  );
}

export default DataNodeItem;
