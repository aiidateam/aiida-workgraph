import * as React from "react";
import { Button, Switch } from "antd";
import { createRender, useModel, useModelState } from "@anywidget/react";
import { createEditor, loadJSON, removeNode, addNode, addLink, removeLink } from "./default_rete";
import "./widget.css";

export function useRete<T extends { destroy(): void }>(
  create: (el: HTMLElement, settings: any) => Promise<T>, settings: any
) {
  const [container, setContainer] = React.useState<null | HTMLElement>(null);
  const editorRef = React.useRef<T>();
  const [editor, setEditor] = React.useState<T | null>(null);
  const [isInitialized, setIsInitialized] = React.useState(false); // Track initialization
  const ref = React.useRef(null);

  React.useEffect(() => {
    if (container) {
      if (editorRef.current) {
        editorRef.current.destroy();
        container.innerHTML = '';
      }
      create(container, settings).then((value) => {
        editorRef.current = value;
        setEditor(value);
        setIsInitialized(true); // Set initialized to true
        window.editor = value;
      });
    }
  }, [container, create]);

  React.useEffect(() => {
    return () => {
      if (editorRef.current) {
        editorRef.current.destroy();
      }
    };
  }, []);

  React.useEffect(() => {
    if (ref.current) {
      setContainer(ref.current);
    }
  }, [ref.current]);

  return [ref, editor, isInitialized] as const;
}


// Function to change the title color based on the state
const changeTitleColor = (ref: any, stateData: any) => {
  if (ref && ref.current) {
    // Find all elements with data-testid="node"
    const nodeElements = ref.current.querySelectorAll('[data-testid="node"]');
    // Iterate through the elements and update title colors based on state
    nodeElements.forEach((nodeElement:any) => {
      const titleElement = nodeElement.querySelector('[data-testid="title"]')  as HTMLElement;
      const nodeName = titleElement.textContent;
      if (nodeName) {
        const nodeState = stateData[nodeName];
        if (nodeState === 'FINISHED') {
          titleElement.style.background = 'green';
        } else if (nodeState === 'RUNNING') {
          titleElement.style.background = 'orange';
        } else if (nodeState === 'CREATED') {
          titleElement.style.background = 'blue';
        } else if (nodeState === 'WAITING') {
          titleElement.style.background = 'purple'; // Change to the desired color for "waiting"
        } else if (nodeState === 'KILLED') {
          titleElement.style.background = 'red'; // Change to the desired color for "killed"
        // } else if (nodeState === 'paused') {
          // titleElement.style.background = 'purple'; // Change to the desired color for "paused"
        } else {
          // Handle any other states or provide a default color
          titleElement.style.background = 'gray'; // Change to the desired default color
        }
      }
    }
  )};
}

const render = createRender(() => {
  const [value, setValue] = useModelState<any>("value");
  const [style, setStyle] = useModelState<any>("style");
  const [settings, setSettings] = useModelState<any>("settings");
  const [states, setStates] = useModelState<any>("states");
  const [positions, setPositions] = useModelState<any>("positions");
  const [ref, editor, isInitialized] = useRete(createEditor, settings);
  const model = useModel();

  function handle_custom_msg(msg: any) {
    console.log(msg.data);
    if (!editor) {
      return;
    }
    let node;
    switch (msg.type) {
      case "test":
        node = editor.editor.getNodes()[0];
        node.label = "hello";
        editor.area.update('node', node.id)
        break;
      case "add_node":
        addNode(editor.editor, editor.area, msg.data)
        // editor.layout(true);
        break;
      case "delete_node":
        removeNode(editor.editor, msg.data.name)
        break;
      case "add_link":
        console.log("Adding link", msg.data);
        addLink(editor.editor, editor.area, editor.layout, msg.data)
        break;
      case "delete_link":
        console.log("Removing link", msg.data);
        removeLink(editor.editor, msg.data)
        break;
      default:
        break;
    }
  }

  React.useEffect(() => {
    if (!editor) {
      return;
    }
    //Add event linstener, the nodetranslated event is fired:
    editor.area.addPipe(context => {
      if (context.type === 'nodetranslated') {
        const view = editor.area.nodeViews.get(context.data.id)
        if (view) {
          console.log("Node translated", view.position)
          const node = editor.editor.getNode(context.data.id);
          console.log("new positions: ", {[node.label]: [view.position.x, view.position.y]})
          setPositions({[node.label]: [view.position.x, view.position.y]})
        }
      }
      return context
    })
  }, [isInitialized]);

  React.useEffect(() => {
    model.on("msg:custom", handle_custom_msg);
    return () => model.off("msg:custom", handle_custom_msg);
  }, [model, isInitialized]);

  React.useEffect(() => {
    console.log("Updating node state", states);
    changeTitleColor(ref, states)
  }, [states, isInitialized]);

  React.useEffect(() => {
    console.log("Value changed", value);
    if (!editor) {
      console.log("Editor not ready");
      return;
    };
    // same uuid skip
    if (editor.editor.uuid && editor.editor.uuid === value.uuid) {
      console.log("Same uuid, skip");
      return;
    }
    loadJSON(editor.editor, editor.area, editor.layout, value);
  }, [isInitialized, value]);

  return (
    <div className="App">
      <div>
        <Button type="primary" onClick={() => editor?.layout(true)}>Arrange</Button>
        {/* <Button type="primary" onClick={() => editor?.layout(true)}>Save</Button> */}
      </div>
      <div ref={ref} className="rete" style={style}></div>
    </div>
  );
});

export default { render };
