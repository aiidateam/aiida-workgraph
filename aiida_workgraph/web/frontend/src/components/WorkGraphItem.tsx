import React, { useState, useEffect, useRef, useMemo} from 'react';
import { useParams } from 'react-router-dom';
import '../App.css';
import '../rete.css';
import { createEditor } from '../rete/default';
import { Button, Switch } from "antd";
import WorkGraphIndicator from './WorkGraphIndicator'; // Import the WorkGraphIndicator component
import WorktreeSummary from './WorkGraphSummary';
import WorkGraphLog from './WorkGraphLog';
import NodeDetails from './NodeDetails';
import NodeDurationGraph from './WorkGraphDuration'
import {
  PageContainer,
  EditorContainer,
  LayoutAction,
  TopMenu,
  EditorWrapper,
} from './WorkGraphItemStyles'; // Import your styles
import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';


// Extend the Window interface
declare global {
  interface Window {
    editor?: any; // You can replace `any` with a more specific type if available
  }
}



/* Modify the useRete function to support passing workgraph data to createEditor */
export function useRete<T extends { destroy(): void }>(
  create: (el: HTMLElement, data: any) => Promise<T>,
  workgraphData: any
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
      create(container, workgraphData).then((value) => {
        editorRef.current = value;
        setEditor(value);
        window.editor = value;
      });
    }
  }, [container, create, workgraphData]); // Add workgraphData as a dependency

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



function WorkGraphGraph() {
  const { pk } = useParams();
  const [workgraphData, setWorktreeData] = useState({ summary: {}, nodes: {}, links: [], pk: [] });
  const [ref, editor] = useRete(createEditor, workgraphData);
  const [selectedNode, setSelectedNode] = useState({ metadata: [], executor: '' });
  const [showNodeDetails, setShowNodeDetails] = useState(false);
  const [workgraphHierarchy, setWorktreeHierarchy] = useState([]);
  const [selectedView, setSelectedView] = useState('Editor');
  const [realtimeSwitch, setRealtimeSwitch] = useState(false); // State to manage the realtime switch

  // Fetch state data from the backend
  const fetchStateData = async () => {
    try {
      const response = await fetch(`http://localhost:8000/api/workgraph-state/${pk}`);
      if (!response.ok) {
        throw new Error('Failed to fetch state data');
      }
      const data = await response.json();
      // Call changeTitleColor here to update title colors based on the new state data
      changeTitleColor(data);
    } catch (error) {
      console.error('Error fetching state data:', error);
    }
  };

  // Function to change the title color based on the state
  const changeTitleColor = (stateData: any) => {
    if (editor && editor.editor) {
      const nodeEditor = editor.editor;
      const nodes = nodeEditor.getNodes();
      // Find all elements with data-testid="node"
      const nodeElements = document.querySelectorAll('[data-testid="node"]');
      // Iterate through the elements and update title colors based on state
      nodeElements.forEach((nodeElement) => {
        const titleElement = nodeElement.querySelector('[data-testid="title"]')  as HTMLElement;
        const nodeName = titleElement.textContent;
        if (nodeName && nodeName in stateData) {
          const nodeState = stateData[nodeName].state;
            if (nodeState === 'FINISHED') {
              titleElement.style.background = 'green';
            } else if (nodeState === 'RUNNING') {
              titleElement.style.background = 'orange';
            } else if (nodeState === 'CREATED') {
              titleElement.style.background = 'blue';
            } else if (nodeState === 'PLANNED') {
              titleElement.style.background = 'gray';
            } else if (nodeState === 'WAITING') {
              titleElement.style.background = 'purple'; // Change to the desired color for "waiting"
            } else if (nodeState === 'KILLED') {
              titleElement.style.background = 'pink'; // Change to the desired color for "killed"
            } else if (nodeState === 'PAUSED') {
              titleElement.style.background = 'yellow'; // Change to the desired color for "paused"
            } else if (nodeState === 'FAILED') {
              titleElement.style.background = 'red'; // Change to the desired color for "failed"
            } else {
              // Handle any other states or provide a default color
              titleElement.style.background = 'lightblue'; // Change to the desired default color
            }
        }
      }
      )};
  }

  // Setup interval for fetching real-time data when the switch is turned on
  useEffect(() => {
    let intervalId: NodeJS.Timeout;
    if (realtimeSwitch) {
      // Fetch data initially
      fetchStateData();
      // Set up an interval to fetch data every 5 seconds
      intervalId = setInterval(fetchStateData, 5000);
    }
    return () => {
      if (intervalId) {
        clearInterval(intervalId); // Clear the interval when the component unmounts or the switch is turned off
      }
    };
  }, [realtimeSwitch, pk, editor]); // Depend on the realtimeSwitch, pk, and editor


  // Fetch workgraph data from the API
  useEffect(() => {
    fetch(`http://localhost:8000/api/workgraph/${pk}`)
      .then((response) => response.json())
      .then((data) => {
        setWorktreeData(data);
        // Set the workgraph hierarchy here based on your data
        setWorktreeHierarchy(data.parent_workgraphs);
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
            const response = await fetch(`http://localhost:8000/api/workgraph/${pk}/${node.label}`);
            if (!response.ok) {
              throw new Error('Failed to fetch data');
            }

            const data = await response.json();
            // Update the component state with the fetched data
            setSelectedNode(data);
          } catch (error) {
            console.error('Error fetching data:', error);
          }

          setShowNodeDetails(true);
        }
        return context;
      };

      editor.area.addPipe(handleNodePick);
      /* Add arrange node, maybe there is a better plance to add */
      editor?.layout(true)


      // Cleanup function to remove the event listener
      // return () => {
      // editor.area.removePipe(handleNodePick);
      // };
    }
  }, [editor]); // Depend on a stable reference of `editor`, `pk`, and `node_name`

  const handleNodeDetailsClose = () => {
    setShowNodeDetails(false);
  };

  // Memoize the editor to prevent re-creation
  const editorComponent = useMemo(() => (
      <div ref={ref} style={{ height: 'calc(100% - 2em)', width: '100%' }}></div>
  ), [workgraphHierarchy, editor, showNodeDetails, selectedNode]); // Specify dependencies


  const handleTaskAction = async (action: string) => {
    if (editor && editor.editor) {
        const selectedNodes = editor.editor.getNodes().filter((node: any) => node.selected);
        const nodeNames = selectedNodes.map((node: any) => node.label);
        console.log(nodeNames); // Good for debugging

        try {
            const response = await fetch(`http://localhost:8000/api/workgraph/tasks/${action}/${pk}`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(nodeNames), // Correct the body to match expected structure
            });

            if (!response.ok) throw new Error(`Failed to perform ${action}: ${response.statusText}`);

            const data = await response.json();
            console.log(data.message); // Display backend response message
            toast.success(`${action} action performed successfully on nodes.`);
        } catch (error: any) {
            console.error('Error performing node action:', error);
            toast.error(`Error performing ${action}: ${error.message}`);
        }
    } else {
        toast.error("No nodes selected or editor is not available");
    }
};

  const handlePause = () => handleTaskAction('pause');
  const handlePlay = () => handleTaskAction('play');
  const handleKill = () => handleTaskAction('kill');

  return (
      <PageContainer>
        <TopMenu>
          <Button onClick={() => setSelectedView('Editor')}>Editor</Button>
          <Button onClick={() => setSelectedView('Summary')}>Summary</Button>
          <Button onClick={() => setSelectedView('Log')}>Log</Button>
          <Button onClick={() => setSelectedView('Time')}>Time</Button>
        </TopMenu>
          {selectedView === 'Summary' && <WorktreeSummary summary={workgraphData.summary} />}
          {selectedView === 'Log' && <WorkGraphLog id={pk} />}
          {selectedView === 'Time' && <NodeDurationGraph id={pk}/>}
          <EditorWrapper visible={selectedView === 'Editor'}>
          <WorkGraphIndicator parentWorktrees={workgraphHierarchy} />
            <EditorContainer>
              <LayoutAction>
              <div>
                <Switch
                  checked={realtimeSwitch}
                  onChange={(checked) => setRealtimeSwitch(checked)}
                  style={{ marginRight: '10px' }}
                />
                <label>Real-time state</label>
              </div>
              <div>
                <ToastContainer />
                <Button onClick={() => editor?.layout(true)}>Arrange</Button>
                <Button onClick={handlePause}>Pause</Button>
                <Button onClick={handlePlay}>Play</Button>
                <Button onClick={handleKill}>Kill</Button>
              </div>
              </LayoutAction>
              {showNodeDetails && (
              <NodeDetails selectedNode={selectedNode} onClose={handleNodeDetailsClose} setShowNodeDetails={setShowNodeDetails} />
            )}
            </EditorContainer>
            {editorComponent}
          </EditorWrapper>
      </PageContainer>
  );
}

export default WorkGraphGraph;
