// WorkGraphLog.js
import styled from "styled-components";


export const WorktreeLogStyle = styled.div`
  .log-section {
  border: 1px solid #ddd;
  padding: 1em;
  overflow-x: auto; /* Add horizontal scrollbar */
  overflow-y: auto; /* Add vertical scrollbar */
  font-family: monospace;
  white-space: pre;
  font-size: 1.2em;
  color: #444;
  line-height: 1.4;
  text-align: left;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
}

.log-content {
  flex-grow: 1;
}
`;

function WorkGraphLog({ logs }) {
  return (
    <WorktreeLogStyle>

    <div className="log-section">
      <h3>Log Information</h3>
      <div className="log-content">
        {logs.map((log, index) => (
          <div key={index}>{log}</div>
        ))}
      </div>
    </div>
    </WorktreeLogStyle>
  );
}

export default WorkGraphLog;
