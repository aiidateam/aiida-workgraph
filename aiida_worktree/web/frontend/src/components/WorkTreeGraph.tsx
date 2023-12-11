import React, { useState, useEffect, useRef } from 'react';
import { useParams } from 'react-router-dom';
import logo from '../logo.svg';
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

const LayoutAction = styled.div`
  position: absolute;
  top: 1em;
  right: 1em;
  display: flex;
  align-items: center;
  gap: 0.5em;
`;

function WorkTreeGraph() {
  const { pk } = useParams(); // Assuming 'pk' is the parameter in the URL
  const [worktreeData, setWorktreeData] = useState({ nodes: {}, links: [] });
  const [animate, setAnimate] = useState(true);

  // Fetch worktree data from the API
  useEffect(() => {
    fetch(`http://localhost:8000/api/worktree/${pk}`)
      .then(response => response.json())
      .then(data => setWorktreeData(data))
      .catch(error => console.error("Error fetching data:", error));
  }, [pk]); // Dependency array, re-fetch if 'pk' changes

  // Only initialize the editor if worktreeData is available
  const [ref, editor] = useRete(createEditor, worktreeData);

  return (
    <div className="App">
      <LayoutAction>
        Animate
        <Switch checked={animate} onChange={setAnimate} />
        <Button onClick={() => editor?.layout(animate)}>Layout</Button>
      </LayoutAction>
      <div ref={ref} style={{ height: "100vh", width: "100vw" }}></div>
    </div>
  );
}

export default WorkTreeGraph;
