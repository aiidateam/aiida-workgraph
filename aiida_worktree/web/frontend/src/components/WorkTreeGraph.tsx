import React, { useState, useEffect, useRef } from 'react';
import { useParams } from 'react-router-dom';
import '../App.css';
import '../rete.css';
import { createEditor } from '../rete/default';
import { Button, Switch } from "antd";
import WorktreeIndicator from './WorktreeIndicator'; // Import the WorktreeIndicator component
import WorktreeSummary from './WorktreeSummary';
import WorkTreeLog from './WorkTreeLog';
import NodeDetails from './NodeDetails';
import {
  PageContainer,
  WorktreeInfo,
  EditorContainer,
  LayoutAction,
} from './WorkTreeGraphStyles'; // Import your styles


/* Modify the useRete function to support passing worktree data to createEditor */
export function useRete<T extends { destroy(): void }>(
  create: (el: HTMLElement, data: any) => Promise<T>,
  worktreeData: any
) {
  const [container, setContainer] = useState<null | HTMLElement>(null);
  const editorRef = useRef<T>();
  const [editor, setEditor] = useState<T | null>(null);
  const ref = useRef(null);

  useEffect(() => {
    if (container) {
      if (editorRef.current) {
        editorRef.current.destroy();
        container.innerHTML = '';
      }
      create(container, worktreeData).then((value) => {
        editorRef.current = value;
        setEditor(value);
      });
    }
  }, [container, create, worktreeData]); // Add worktreeData as a dependency

  useEffect(() => {
    return () => {
      if (editorRef.current) {
        editorRef.current.destroy();
      }
    };
  }, []);

  useEffect(() => {
    if (ref.current) {
      setContainer(ref.current);
    }
  }, [ref.current]);

  return [ref, editor] as const;
}




function WorkTreeGraph() {
  const { pk } = useParams();
  const [worktreeData, setWorktreeData] = useState({ summary: [], nodes: {}, links: [], logs: [] });
  const [animate, setAnimate] = useState(true);
  const [ref, editor] = useRete(createEditor, worktreeData);
  const [selectedNode, setSelectedNode] = useState({ metadata: [], executor: '' });
  const [showNodeDetails, setShowNodeDetails] = useState(false);
  const [worktreeHierarchy, setWorktreeHierarchy] = useState([]);

  // Fetch worktree data from the API
  useEffect(() => {
    fetch(`http://localhost:8000/api/worktree/${pk}`)
      .then((response) => response.json())
      .then((data) => {
        setWorktreeData(data);
        // Set the worktree hierarchy here based on your data
        console.log(data.parent_worktrees);
        setWorktreeHierarchy(data.parent_worktrees);
      })
      .catch((error) => console.error('Error fetching data:', error));
  }, [pk]); // Only re-run when `pk` changes

  // Setup editor event listener
  useEffect(() => {
    if (editor) {
      const handleNodePick = async (context: any) => {
        if (!context || typeof context !== 'object' || !('type' in context)) return context;

        if (context.type === 'nodepicked') {
          const pickedId = context.data.id;
          const node = editor.editor.getNode(pickedId);

          try {
            // Fetch data from the backend
            const response = await fetch(`http://localhost:8000/api/worktree/${pk}/${node.label}`);
            if (!response.ok) {
              throw new Error('Failed to fetch data');
            }

            const data = await response.json();
            // Assuming data contains the details you need from the backend
            // Update your component state with the fetched data
            setSelectedNode(data);
          } catch (error) {
            console.error('Error fetching data:', error);
          }

          setShowNodeDetails(true);
        }
        return context;
      };

      editor.area.addPipe(handleNodePick);

      // Cleanup function to remove the event listener
      // return () => {
      // editor.area.removePipe(handleNodePick);
      // };
    }
  }, [editor]); // Depend on a stable reference of `editor`, `pk`, and `node_name`

  const handleNodeDetailsClose = () => {
    setShowNodeDetails(false);
  };

  return (
    <div className="App">
      <PageContainer>
        <WorktreeInfo>
          <WorktreeSummary summary={worktreeData.summary} />
          <WorkTreeLog logs={worktreeData.logs} />
        </WorktreeInfo>
        <EditorContainer>
          <WorktreeIndicator parentWorktrees={worktreeHierarchy} />
          <LayoutAction>
            <label>
              Animate
              <Switch checked={animate} onChange={setAnimate} />
            </label>
            <Button onClick={() => editor?.layout(animate)}>Arrange</Button>
          </LayoutAction>
          <div ref={ref} style={{ height: 'calc(100% - 2em)', width: '100%' }}></div>
          {showNodeDetails && (
            <NodeDetails selectedNode={selectedNode} onClose={handleNodeDetailsClose} setShowNodeDetails={setShowNodeDetails} />
          )}
        </EditorContainer>
      </PageContainer>
    </div>
  );
}

export default WorkTreeGraph;
