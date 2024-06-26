// WorkGraphLog.js
import styled from "styled-components";
import { useEffect, useState } from "react";


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

function WorkGraphLog({ id }) {
  const [fetchedLogs, setFetchedLogs] = useState([]);

  useEffect(() => {
    fetchLogs(); // Fetch logs immediately

    const interval = setInterval(() => {
      fetchLogs();
    }, 4000);

    return () => {
      clearInterval(interval);
    };
  }, []);

  const fetchLogs = async () => {
    try {
      const response = await fetch(`http://localhost:8000/api/workgraph-logs/${id}`);
      const data = await response.json();
      setFetchedLogs(data);
    } catch (error) {
      console.error("Error fetching logs:", error);
    }
  };

  return (
    <WorktreeLogStyle>
      <div className="log-section">
        <h3>Log Information</h3>
        <div className="log-content">
          {fetchedLogs.map((log, index) => (
            <div key={index}>{log}</div>
          ))}
        </div>
      </div>
    </WorktreeLogStyle>
  );
}

export default WorkGraphLog;
