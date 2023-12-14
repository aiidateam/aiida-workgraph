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
  export const InputsCode = styled(SyntaxHighlighter)`
  width: 100%;
  max-width: 100%;
  max-height: 300px;
  overflow-x: auto;
  white-space: pre;
  margin-top: 10px;
  border: 1px solid #ccc;
  border-radius: 4px;
  padding: 10px;
  background-color: #f7f7f7;
  font-family: monospace;
`;

function WorktreeSummary({ summary }) {
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
        <strong>Inputs:</strong>
      </div>
      <InputsCode language="python" style={dark}>
        {/* Display args and inputs here */}
        {summary.inputs}
      </InputsCode>
      <div>
        <strong>Outputs:</strong>
      </div>
      <InputsCode language="python" style={dark}>
        {/* Display args and inputs here */}
        {summary.outputs}
      </InputsCode>
    </div>
    </WorktreeInfoStyle>
  );
}

export default WorktreeSummary;
