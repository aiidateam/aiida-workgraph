import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Home from './components/Home';
import WorkTreeTable from './components/WorkTreeTable';
import Node from './components/Node';
import WorkTreeItem from './components/WorkTreeItem';
import Settings from './components/Settings';
import Layout from './components/Layout'; // Import the Layout component
import './App.css';

function App() {
  return (
    <Router>
      <div className="App">
        <Layout> {/* Wrap the routes with the Layout component */}
          <Routes>
            <Route path="/worktree" element={<WorkTreeTable />} />
            <Route path="/node" element={<Node />} />
            <Route path="/settings" element={<Settings />} />
            <Route path="/" element={<Home />} />
            <Route path="/worktree/:pk" element={<WorkTreeItem />} />
          </Routes>
        </Layout>
      </div>
    </Router>
  );
}

export default App;
