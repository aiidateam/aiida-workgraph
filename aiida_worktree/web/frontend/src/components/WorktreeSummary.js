// WorktreeSummary.js
import React from 'react';

function WorktreeSummary({ summary }) {
  return (
    <div>
      <h2>Summary</h2>
      <div className="info-table">
        {summary.map(([property, value]) => (
          <div className="info-row" key={property}>
            <div className="property">{property}</div>
            <div className="value">{value}</div>
          </div>
        ))}
      </div>
    </div>
  );
}

export default WorktreeSummary;
