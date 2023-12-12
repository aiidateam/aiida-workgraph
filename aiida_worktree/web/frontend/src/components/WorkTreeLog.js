// WorkTreeLog.js
import React from 'react';

function WorkTreeLog({ logs }) {
  return (
    <div className="log-section">
      <h3>Log Information</h3>
      {logs.map((log, index) => (
        <div key={index}>{log}</div>
      ))}
    </div>
  );
}

export default WorkTreeLog;
