// WorktreeSummary.js
import styled from "styled-components";

export const WorktreeInfoStyle = styled.div`
  padding: 1em;
  overflow-y: auto;
  box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
  box-sizing: border-box;
  display: flex;
  flex-direction: column;
  height: 100%;

  h2 {
    margin-bottom: 0.5em;
    color: #333; // Darker color for headers
    font-size: 1.2em;
  }

  .info-table {
    flex-grow: 1; // Allow this section to take available space
    overflow-y: auto; // Make only this section scrollable if needed
    margin-bottom: 1em;

    .info-row {
      display: flex;
      border-bottom: 1px solid #eee;
      padding: 0.5em 0;

      .property {
        width: 40%; // Adjust for better alignment
        font-weight: bold;
        text-align: left;
        font-size: 1.2em;
        color: #555; // Slightly darker for better readability
      }

      .value {
        width: 60%; // Adjust accordingly
        text-align: left;
        font-size: 1.2em;
        color: #666; // Slightly lighter to differentiate from property
      }
    }
  }
  `;


function WorktreeSummary({ summary }) {
  return (
    <WorktreeInfoStyle>
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
    </WorktreeInfoStyle>
  );
}

export default WorktreeSummary;
