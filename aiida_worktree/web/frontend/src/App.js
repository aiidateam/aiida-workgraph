import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Home from './components/Home';
import WorkTreeTable from './components/WorkTreeTable';
import DataNodeTable from './components/DataNodeTable';
import WorkTreeItem from './components/WorkTreeItem';
import DataNodeItem from './components/DataNodeItem';
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
            <Route path="/datanode" element={<DataNodeTable />} />
            <Route path="/settings" element={<Settings />} />
            <Route path="/" element={<Home />} />
            <Route path="/worktree/:pk" element={<WorkTreeItem />} />
            <Route path="/datanode/:pk" element={<DataNodeItem />} />
          </Routes>
        </Layout>
      </div>
    </Router>
  );
}

export default App;
