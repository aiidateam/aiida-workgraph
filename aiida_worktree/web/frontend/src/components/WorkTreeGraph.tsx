import React, { useState, useEffect, useRef } from 'react';
import { useParams } from 'react-router-dom';
import '../App.css';
import '../rete.css';
import { createEditor } from '../rete/default';
import styled from "styled-components";
import { Button, Switch } from "antd";
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { dark } from 'react-syntax-highlighter/dist/esm/styles/prism'; // Correct import for 'dark' style

/* Modify the useRete function to support passing worktree data to createEditor */
export function useRete<T extends { destroy(): void }>(
  create: (el: HTMLElement, data: any) => Promise<T>,
  worktreeData: any
) {
  const [container, setContainer] = useState<null | HTMLElement>(null);
  const editorRef = useRef<T>();
  const [editor, setEditor] = useState<T | null>(null);
  const ref = useRef(null);

  useEffect(() => {
    if (container) {
      if (editorRef.current) {
        editorRef.current.destroy();
        container.innerHTML = '';
      }
      create(container, worktreeData).then((value) => {
        editorRef.current = value;
        setEditor(value);
      });
    }
  }, [container, create, worktreeData]); // Add worktreeData as a dependency

  useEffect(() => {
    return () => {
      if (editorRef.current) {
        editorRef.current.destroy();
      }
    };
  }, []);

  useEffect(() => {
    if (ref.current) {
      setContainer(ref.current);
    }
  }, [ref.current]);

  return [ref, editor] as const;
}


const PageContainer = styled.div`
  display: flex;
  height: 100vh;
  width: 100vw;
  background-color: #f4f4f4;
`;

const WorktreeInfo = styled.div`
  width: 20%;
  padding: 1em;
  overflow-y: auto;
  background-color: #fff;
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
        font-size: 0.9em;
        color: #555; // Slightly darker for better readability
      }

      .value {
        width: 60%; // Adjust accordingly
        text-align: left;
        font-size: 0.8em;
        color: #666; // Slightly lighter to differentiate from property
      }
    }
  }

  .log-section {
    background-color: #f9f9f9;
    border: 1px solid #ddd;
    padding: 1em;
    overflow-y: auto;
    max-height: 200px;
    font-family: monospace;
    white-space: pre-wrap;
    font-size: 0.9em;
    color: #444;
    line-height: 1.4;
    text-align: left; // Ensure text is left-aligned
    display: flex;
    flex-direction: column;
    align-items: flex-start; // Aligns children (log entries) to the start (left)
  }
`;


const EditorContainer = styled.div`
  width: 80%;
  padding: 1em;
  box-sizing: border-box;
`;

const LayoutAction = styled.div`
  position: absolute;
  top: 1em;
  right: 1em;
  display: flex;
  align-items: center;
  gap: 0.5em;
  z-index: 10;
`;


const NodeDetailsPanel = styled.div`
  position: absolute;
  top: 0;
  right: 0;
  background-color: #fff;
  box-shadow: -2px 0 5px rgba(0, 0, 0, 0.1);
  width: 300px;
  height: 100%;
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
  flex-grow: 1; /* Allow this section to take available space */
  overflow-y: auto; /* Make only this section scrollable if needed */
  margin-bottom: 1em;
`;

const NodeDetailRow = styled.div`
  display: flex;
  border-bottom: 1px solid #eee;
  padding: 0.5em 0;
`;

const NodeDetailProperty = styled.div`
  width: 40%; /* Adjust for better alignment */
  font-weight: bold;
  text-align: left;
  font-size: 0.9em;
  color: #555; /* Slightly darker for better readability */
`;

const NodeDetailValue = styled.div`
  width: 60%; /* Adjust accordingly */
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


function WorkTreeGraph() {
  const { pk } = useParams();
  const [worktreeData, setWorktreeData] = useState({ summary: [], nodes: {}, links: [] , logs: []});
  const [animate, setAnimate] = useState(true);
  const [ref, editor] = useRete(createEditor, worktreeData);
  const [selectedNode, setSelectedNode] = useState({metadata: [], executor: ""});
  const [showNodeDetails, setShowNodeDetails] = useState(false);


  // Fetch worktree data from the API
  useEffect(() => {
    fetch(`http://localhost:8000/api/worktree/${pk}`)
      .then(response => response.json())
      .then(data => setWorktreeData(data))
      .catch(error => console.error("Error fetching data:", error));
  }, [pk]); // Only re-run when `pk` changes

  // Setup editor event listener
  useEffect(() => {
    if (editor) {
      const handleNodePick = async (context: any) => {
        if (!context || typeof context !== 'object' || !('type' in context)) return context;

        if (context.type === 'nodepicked') {
          const pickedId = context.data.id;
          const node = editor.editor.getNode(pickedId);

          try {
            // Fetch data from the backend
            const response = await fetch(`http://localhost:8000/api/worktree/${pk}/${node.label}`);
            if (!response.ok) {
              throw new Error('Failed to fetch data');
            }

            const data = await response.json();
            // Assuming data contains the details you need from the backend
            // Update your component state with the fetched data
            setSelectedNode(data);
          } catch (error) {
            console.error('Error fetching data:', error);
          }

          setShowNodeDetails(true);
        }
        return context;
      };

      editor.area.addPipe(handleNodePick);

      // Cleanup function to remove the event listener
      // return () => {
        // editor.area.removePipe(handleNodePick);
      // };
    }
  }, [editor]); // Depend on a stable reference of `editor`, `pk`, and `node_name`


  return (
    <div className="App">
      <PageContainer>
        <WorktreeInfo>
          <h2>Summary</h2>
          <div className="info-table">
            {/* Iterate over worktreeData.summary */}
            {worktreeData.summary.map(([property, value]) => (
              <div className="info-row" key={property}>
                <div className="property">{property}</div>
                <div className="value">{value}</div>
              </div>
            ))}
          </div>

          <div className="log-section">
            <h3>Log Information</h3>
            {/* Display each log entry */}
            {worktreeData.logs.map((log, index) => (
              <div key={index}>{log}</div>
            ))}
          </div>
        </WorktreeInfo>
        <EditorContainer>
          <LayoutAction>
            <label>
              Animate
              <Switch checked={animate} onChange={setAnimate} />
            </label>
            <Button onClick={() => editor?.layout(animate)}>Layout</Button>
          </LayoutAction>
          <div ref={ref} style={{ height: "calc(100% - 2em)", width: "100%" }}></div>
          {/* Node Details Panel */}
          {showNodeDetails && (
            <NodeDetailsPanel>
            <CloseButton onClick={() => setShowNodeDetails(false)}>Close</CloseButton>
            <NodeDetailsTitle>Node Details</NodeDetailsTitle>
            {/* Display the selected node's details */}
            {selectedNode && (
              <NodeDetailsTable>
                {/* Render the details of the selected node */}
                {selectedNode.metadata.map(([property, value]) => (
                  <NodeDetailRow key={property}>
                    <NodeDetailProperty>{property}</NodeDetailProperty>
                    <NodeDetailValue>{value}</NodeDetailValue>
                  </NodeDetailRow>
                ))}
              </NodeDetailsTable>
            )}
            <div>
              <strong>Executor:</strong>
            </div>
            {/* Highlight and format the Python code */}
            <PythonCode language="python" style={dark}>
              {selectedNode.executor}
            </PythonCode>
            {/* Add more details as needed */}
          </NodeDetailsPanel>
          )}
        </EditorContainer>
      </PageContainer>
    </div>
  );
}

export default WorkTreeGraph;
