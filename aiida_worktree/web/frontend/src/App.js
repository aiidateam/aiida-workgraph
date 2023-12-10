import React from 'react';
import { BrowserRouter as Router, Route, Link, Routes } from 'react-router-dom';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faHome, faTree, faDotCircle, faCog } from '@fortawesome/free-solid-svg-icons';
import Home from './components/Home';
import WorkTree from './components/WorkTree';
import Node from './components/Node';
import WorkTreeGraph from './components/WorkTreeGraph'; // Import the component for the detail page
import Settings from './components/Settings'; // Import your Settings component
import './App.css';

function App() {
  return (
    <Router>
      <div className="App">
        <div className="sidebar">
          <nav>
            <ul>
              <li><Link to="/"><FontAwesomeIcon icon={faHome} /><span>Home</span></Link></li>
              <li><Link to="/worktree"><FontAwesomeIcon icon={faTree} /><span>WorkTree</span></Link></li>
              <li><Link to="/node"><FontAwesomeIcon icon={faDotCircle} /><span>Node</span></Link></li>
              <li><Link to="/settings"><FontAwesomeIcon icon={faCog} /><span>Settings</span></Link></li>
            </ul>
          </nav>
        </div>
        <div className="content">
          <Routes>
            <Route path="/worktree" element={<WorkTree />} />
            <Route path="/node" element={<Node />} />
            <Route path="/settings" element={<Settings />} />
            <Route path="/" element={<Home />} />
            <Route path="/worktree/:pk" element={<WorkTreeGraph />} />
          </Routes>
        </div>
      </div>
    </Router>
  );
}

export default App;
