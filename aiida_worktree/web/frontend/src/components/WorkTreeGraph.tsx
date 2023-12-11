import React, { useState, useEffect, useRef } from 'react';
import { useParams } from 'react-router-dom';
import '../App.css';
import '../rete.css';
import { createEditor } from '../rete/default';
import styled from "styled-components";
import { Button, Switch } from "antd";

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


function WorkTreeGraph() {
  const { pk } = useParams();
  const [worktreeData, setWorktreeData] = useState({ summary: [], nodes: {}, links: [] , logs: []});
  const [animate, setAnimate] = useState(true);
  const [ref, editor] = useRete(createEditor, worktreeData);

  // Fetch worktree data from the API
  useEffect(() => {
    fetch(`http://localhost:8000/api/worktree/${pk}`)
      .then(response => response.json())
      .then(data => {
        // Assuming 'data' contains the 'summary' field
        setWorktreeData(data);
      })
      .catch(error => console.error("Error fetching data:", error));
  }, [pk]);

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
        </EditorContainer>
      </PageContainer>
    </div>
  );
}

export default WorkTreeGraph;
