import React from 'react';
import styled from 'styled-components';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { dark } from 'react-syntax-highlighter/dist/esm/styles/prism'; // Correct import for 'dark' style
import { useNavigate } from 'react-router-dom'; // Use the useNavigate hook

const WorktreeButton = styled.button`
  padding: 10px;
  background-color: #007bff;
  color: #fff;
  border: none;
  border-radius: 4px;
  cursor: pointer;
  margin-top: 10px;

  &:disabled {
    background-color: #ccc; /* Gray color for disabled state */
    cursor: not-allowed; /* Change cursor to indicate it's not clickable */
    color: #666; /* Optional: change text color for disabled state */
  }
`;

const NodeDetailsPanel = styled.div`
  position: absolute;
  top: 0;
  right: 0;
  background-color: #fff;
  box-shadow: -2px 0 5px rgba(0, 0, 0, 0.1);
  width: 25%;
  height: 100vh;
  padding: 20px;
  box-sizing: border-box;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  overflow-y: auto;
  z-index: 20;
  border-left: 1px solid #ddd; /* Add a border to separate it from the editor */
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

const NodeDetailRow = styled.div`
  display: flex;
  border-bottom: 1px solid #eee;
  padding: 0.5em 0;
`;

const NodeDetailProperty = styled.div`
  width: 50%; /* Adjust for better alignment */
  font-weight: bold;
  text-align: left;
  font-size: 0.9em;
  color: #555; /* Slightly darker for better readability */
`;

const NodeDetailValue = styled.div`
  width: 50%; /* Adjust accordingly */
  text-align: left;
  font-size: 0.8em;
  color: #666; /* Slightly lighter to differentiate from property */
`;

const CloseButton = styled.button`
  align-self: flex-end;
  margin-bottom: 10px;
  // Style your button as needed
`;

const PythonCode = styled(SyntaxHighlighter)`
  width: 100%;
  max-width: 100%; /* Limit the maximum width to prevent stretching */
  max-height: 300px; /* Limit the maximum height to add a vertical scrollbar */
  overflow-x: auto; /* Add a horizontal scrollbar if needed */
  white-space: pre; /* Preserve whitespace */
  margin-top: 10px;
  border: 1px solid #ccc;
  border-radius: 4px;
  padding: 10px;
  background-color: #f7f7f7;
  font-family: monospace; /* Use a monospace font */
`;

const InputsCode = styled(SyntaxHighlighter)`
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

const StyledLink = styled.span`
  cursor: pointer;
  text-decoration: underline;
  color: blue;
`;

function NodeDetails({ selectedNode, onClose, setShowNodeDetails }) {
  const navigate = useNavigate();

  const handleClose = () => {
    setShowNodeDetails(false);
    onClose();
  };

  const handleWorktreeClick = () => {
    if (selectedNode.node_type === 'graph_builder' && selectedNode.process.pk) {
      navigate(`/workgraph/${selectedNode.process.pk}`);
    }
  };


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

  // Determine if the button should be disabled
  const isButtonDisabled = !selectedNode.process || !selectedNode.process.pk;

  return (
    <NodeDetailsPanel>
      <CloseButton onClick={handleClose}>Close</CloseButton>
      <NodeDetailsTitle>Node Details</NodeDetailsTitle>
      {selectedNode.node_type === 'graph_builder' && (
      <WorktreeButton onClick={handleWorktreeClick} disabled={isButtonDisabled}>
        Go to Worktree
      </WorktreeButton>
      )}
      {selectedNode && (
        <NodeDetailsTable>
          {selectedNode.metadata.map(([property, value]) => (
            <NodeDetailRow key={property}>
              <NodeDetailProperty>{property}</NodeDetailProperty>
              <NodeDetailValue>{value}</NodeDetailValue>
            </NodeDetailRow>
          ))}
        </NodeDetailsTable>
      )}
      <div>
        <NodeDetailsTitle>Inputs:</NodeDetailsTitle>
      </div>
      <NodeDetailsTable>
        <ul style={{ margin: 10, padding: 5, textAlign: 'left' }}>
          {renderInputs(selectedNode.inputs)}
        </ul>
      </NodeDetailsTable>
      <div>
        <NodeDetailsTitle>Outputs:</NodeDetailsTitle>
      </div>
      <NodeDetailsTable>
        <ul style={{ margin: 10, padding: 5, textAlign: 'left' }}>
          {renderInputs(selectedNode.outputs)}
        </ul>
      </NodeDetailsTable>
      <div>
        <NodeDetailsTitle>Executor:</NodeDetailsTitle>
      </div>
      <PythonCode language="python" style={dark}>
        {selectedNode.executor}
      </PythonCode>
    </NodeDetailsPanel>
  );
}

export default NodeDetails;
