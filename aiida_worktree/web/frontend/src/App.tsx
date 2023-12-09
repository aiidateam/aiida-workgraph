import React from 'react';
import { useRete } from 'rete-react-plugin';
import logo from './logo.svg';
import './App.css';
import './rete.css';
import { createEditor } from './rete';

function App() {
  const [ref] = useRete(createEditor)

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
        <div ref={ref} className="rete"></div>
      </header>
    </div>
  );
}

export default App
