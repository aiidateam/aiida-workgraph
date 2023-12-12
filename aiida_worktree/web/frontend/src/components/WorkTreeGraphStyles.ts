import styled from "styled-components";

export const PageContainer = styled.div`
  display: flex;
  flex-direction: column;
  height: 100vh;
  width: 100vw;
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
  top: 1em;
  left: 1em;
  display: flex;
  align-items: center;
  gap: 0.5em;
  z-index: 10;
`;
