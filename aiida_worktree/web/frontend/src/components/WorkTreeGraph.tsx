import React, { useState, useEffect, useRef } from 'react';
import { useParams } from 'react-router-dom';
import logo from '../logo.svg';
import '../App.css';
import '../rete.css';
import { createEditor } from '../rete/default';

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


function WorkTreeGraph() {
  let { pk } = useParams();
  const [worktreeData, setWorktreeData] = useState(null);
  const reteContainerRef = useRef(null);

  useEffect(() => {
    fetch(`http://localhost:8000/worktree/${pk}`)
      .then(response => response.json())
      .then(data => setWorktreeData(data))
      .catch(error => console.error('Error fetching data:', error));
  }, [pk]);

  useEffect(() => {
    if (worktreeData && reteContainerRef.current) {
      createEditor(reteContainerRef.current, worktreeData);
    }
  }, [worktreeData]); // Run this effect when worktreeData changes

  return (
    <div className="App">
      <header className="App-header">
        <img src={logo} className="App-logo" alt="logo" style={{ animation: 'none' }} />
        <a
          className="App-link"
          href="https://rete.js.org"
          target="_blank"
          rel="noopener noreferrer"
        >
          Learn Rete.js
        </a>
        <div ref={reteContainerRef} className="rete"></div>
      </header>
    </div>
  );
}

export default WorkTreeGraph;
