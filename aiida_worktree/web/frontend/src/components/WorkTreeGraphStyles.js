import styled from "styled-components";

export const PageContainer = styled.div`
  display: flex;
  height: 100vh;
  width: 100vw;
  background-color: #f4f4f4;
`;

export const WorktreeInfo = styled.div`
  width: 20%;
  padding: 1em;
  overflow-y: auto;
  background-color: #fff;
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
        font-size: 0.9em;
        color: #555; // Slightly darker for better readability
      }

      .value {
        width: 60%; // Adjust accordingly
        text-align: left;
        font-size: 0.8em;
        color: #666; // Slightly lighter to differentiate from property
      }
    }
  }

  .log-section {
    background-color: #f9f9f9;
    border: 1px solid #ddd;
    padding: 1em;
    overflow-y: auto;
    max-height: 200px;
    font-family: monospace;
    white-space: pre-wrap;
    font-size: 0.9em;
    color: #444;
    line-height: 1.4;
    text-align: left; // Ensure text is left-aligned
    display: flex;
    flex-direction: column;
    align-items: flex-start; // Aligns children (log entries) to the start (left)
  }
`;


export const EditorContainer = styled.div`
  width: 80%;
  padding: 1em;
  box-sizing: border-box;
`;

export const LayoutAction = styled.div`
  position: absolute;
  top: 1em;
  right: 1em;
  display: flex;
  align-items: center;
  gap: 0.5em;
  z-index: 10;
`;
