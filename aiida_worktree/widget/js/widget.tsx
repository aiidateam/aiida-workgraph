import * as React from "react";
import { createRender, useModel, useModelState } from "@anywidget/react";
import { createEditor, createDynamicNode } from "./default_rete";
import "./widget.css";

export function useRete<T extends { destroy(): void }>(
  create: (el: HTMLElement, data: any) => Promise<T>,
  worktreeData: any
) {
  const [container, setContainer] = React.useState<null | HTMLElement>(null);
  const editorRef = React.useRef<T>(); // Ref for storing the actual editor instance
  const [editor, setEditor] = React.useState<T | null>(null);
  const ref = React.useRef(null);
  const editorInstanceRef = React.useRef<T | null>(null); // Ref to hold the current editor for access in callbacks

  React.useEffect(() => {
    if (container) {
      if (editorRef.current) {
        editorRef.current.destroy();
        container.innerHTML = '';
      }
      create(container, worktreeData).then((value) => {
        editorRef.current = value;
        setEditor(value);
        // Expose the editor instance to the window object for debugging
        window.editor = value;
      });
    }
  }, [container, create, worktreeData]);

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

  // Update the ref whenever the editor changes
  React.useEffect(() => {
    editorInstanceRef.current = editor;
  }, [editor]);

  return [ref, editor] as const;
}

const render = createRender(() => {
  const [value, setValue] = useModelState<number>("value");
  const [ref, editor] = useRete(createEditor, value);
  const model = useModel();

  React.useEffect(() => {
    function handle_custom_msg(msg: any) {
      console.log(msg.data);
      // Access the editor through the ref
      if (!editor) {
        return;
      }
      /*
        Message format:
        {
          type: "custom",
          data: {
            node_id: string,
            label: string
          }
        }
        type: could be the following:
        - "test"
        - "add_node"
        - "delete_node"
        - "add_link"
        - "remove_link"
      */
		    // console.log(editor.editor.getNodes()[0].label);
        let node;
        switch (msg.type) {
          case "test":
            node = editor.editor.getNodes()[0];
  		      node.label = "hello";
            editor.area.update('node', node.id)
            break;
          case "add_node":
            node = createDynamicNode(msg.data);
            editor.editor.addNode(node);
            break;
          case "delete_node":
            console.log("delete node: ")
            console.log(editor.editor.nodeMap[msg.data.name]);
            editor.editor.removeNode(editor.editor.nodeMap[msg.data.name].id).then(() => {
              delete editor.editor.nodeMap[msg.data.name];
              console.log("Deleted node successfully");
              editor.area.update();
            });
            break;
          case "add_link":
            const fromNode = editor.editor.editor.nodeMap[link.from_node];
            const toNode = editor.editor.editor.nodeMap[link.to_node];
            if (fromNode && toNode) {
                editor.addConnection(new Connection(fromNode, link.from_socket, toNode, link.to_socket));
            }
            break;
          case "remove_link":
            // Do something
            break;
          default:
            break;
        }
    }
    model.on("msg:custom", handle_custom_msg);
    return () => model.off("msg:custom", handle_custom_msg);
  }, [model, editor]);

  return (
    <div className="App">
      <div ref={ref} className="rete"></div>
    </div>
  );
});

export default { render };
