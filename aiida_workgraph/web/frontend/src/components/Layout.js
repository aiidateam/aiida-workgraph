// src/components/Layout.js
import React from 'react';
import { Link } from 'react-router-dom';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faHome, faTree, faDotCircle, faCog } from '@fortawesome/free-solid-svg-icons';
import './Layout.css'; // Import layout-specific styles

const Layout = ({ children }) => {
  return (
    <div className="App">
      <div className="sidebar">
        <nav>
          <ul>
            <li><Link to="/"><FontAwesomeIcon icon={faHome} /><span>Home</span></Link></li>
            <li><Link to="/workgraph"><FontAwesomeIcon icon={faTree} /><span>WorkGraph</span></Link></li>
            <li><Link to="/datanode"><FontAwesomeIcon icon={faDotCircle} /><span>DataNode</span></Link></li>
            <li><Link to="/settings"><FontAwesomeIcon icon={faCog} /><span>Settings</span></Link></li>
          </ul>
        </nav>
      </div>
      <div className="content">
        {children} {/* This is where your page-specific content will be rendered */}
      </div>
    </div>
  );
};

export default Layout;
