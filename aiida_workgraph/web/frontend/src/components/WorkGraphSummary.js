// WorktreeSummary.js
import styled from "styled-components";
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { dark } from 'react-syntax-highlighter/dist/esm/styles/prism'; // Correct import for 'dark' style

export const WorktreeInfoStyle = styled.div`
  width: 50%;
  padding: 1em;
  overflow-y: auto;
  box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
  box-sizing: border-box;
  display: flex;
  flex-direction: column;
  height: 100%;

  h2 {
    margin-bottom: 0.5em;
    color: #333; // Darker color for headers
    font-size: 1.2em;
  }

  .info-table {
    flex-grow: 1; // Allow this section to take available space
    overflow-y: auto; // Make only this section scrollable if needed
    margin-bottom: 1em;

    .info-row {
      display: flex;
      border-bottom: 1px solid #eee;
      padding: 0.5em 0;

      .property {
        width: 40%; // Adjust for better alignment
        font-weight: bold;
        text-align: left;
        font-size: 1.2em;
        color: #555; // Slightly darker for better readability
      }

      .value {
        width: 60%; // Adjust accordingly
        text-align: left;
        font-size: 1.2em;
        color: #666; // Slightly lighter to differentiate from property
      }
    }
  }
  `;

  const NodeDetailsTitle = styled.h3`
  font-size: 1.2em;
  margin-bottom: 0.5em;
  color: #333; /* Darker color for headers */
`;


  const NodeDetailsTable = styled.div`
  width: 100%;
  flex-grow: 1; /* Allow this section to take available space */
  overflow-y: auto; /* Make only this section scrollable if needed */
  margin-bottom: 1em;
  background-color: #f7f7f7; /* Light gray background for better readability */
`;

function WorktreeSummary({ summary }) {

  const renderInputs = (inputs, depth = 0) => {
    return Object.entries(inputs).map(([key, value]) => {
      const nodeId = Array.isArray(value) ? value[0] : value;
      const nodeType = Array.isArray(value) ? value[1] : null;

      if (Array.isArray(value)) {
        return (
          <li key={key}>
            <span >
              {key}: <a href={`/datanode/${nodeId}`}>{nodeId}</a>
            </span>
          </li>
        );
      } else if (typeof value === 'object') {
        return (
          <li key={key}>
            <span >
              {key}:
            </span>
            <ul>{renderInputs(value, depth + 1)}</ul>
          </li>
        );
      } else {
        return null; // or handle other types if needed
      }
    });
  };

  return (
    <WorktreeInfoStyle>
    <div>
      <h2>Summary</h2>
      <div className="info-table">
        {summary.table.map(([property, value]) => (
          <div className="info-row" key={property}>
            <div className="property">{property}</div>
            <div className="value">{value}</div>
          </div>
        ))}
      </div>
      <div>
        <NodeDetailsTitle>Inputs:</NodeDetailsTitle>
      </div>
      <NodeDetailsTable>
        <ul style={{ margin: 10, padding: 5, textAlign: 'left' }}>
          {renderInputs(summary.inputs)}
        </ul>
      </NodeDetailsTable>
      <div>
        <NodeDetailsTitle>Outputs:</NodeDetailsTitle>
      </div>
      <NodeDetailsTable>
        <ul style={{ margin: 10, padding: 5, textAlign: 'left' }}>
          {renderInputs(summary.outputs)}
        </ul>
      </NodeDetailsTable>
    </div>
    </WorktreeInfoStyle>
  );
}

export default WorktreeSummary;
