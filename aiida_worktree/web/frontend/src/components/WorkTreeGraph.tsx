import React, { useState, useEffect, useRef } from 'react';
import { useParams } from 'react-router-dom';
import { useRete } from 'rete-react-plugin';
import logo from '../logo.svg';
import '../App.css';
import '../rete.css';
import { createEditor } from '../rete/default';

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
          href="https://github.com/superstar54/aiida-worktree"
          target="_blank"
          rel="noopener noreferrer"
        >
        Learn AiiDA-WorkTree
        </a>
        <div ref={reteContainerRef} className="rete"></div>
      </header>
    </div>
  );
}

export default WorkTreeGraph;
