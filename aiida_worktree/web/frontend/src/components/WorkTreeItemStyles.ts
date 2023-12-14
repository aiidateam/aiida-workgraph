import styled from "styled-components";

export const PageContainer = styled.div`
  display: flex;
  flex-direction: column;
  height: 95vh;
  width: 95vw;
  background-color: #f4f4f4;
`;

export const TopMenu = styled.div`
  display: flex;
  justify-content: left; // Center the buttons horizontally
  padding: 1em;
`;

interface EditorWrapperProps {
  visible: boolean;
}

export const EditorWrapper = styled.div<EditorWrapperProps>`
display: ${props => props.visible ? 'block' : 'none'};
flex-grow: 1;
padding: 1em;
box-sizing: border-box;
width: 100%; // Ensure full width
height: 100%; // Ensure full height
`;

export const EditorContainer = styled.div`
  width: 100%;
  padding: 1em;
  box-sizing: border-box;
  flex-grow: 1; // Allows the editor to take up the remaining vertical space
`;

export const LayoutAction = styled.div`
  position: relative;
  display: flex;
  flex-direction: column; /* Removed single quotes */
  align-items: flex-start; /* Removed single quotes */
  top: 1em;
  left: 1em;
  align-items: left;
  gap: 0.5em;
  z-index: 10;
`;
