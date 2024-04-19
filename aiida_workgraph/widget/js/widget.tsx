import * as React from "react";
import { createRender, useModel, useModelState } from "@anywidget/react";
import { createEditor, removeNode, addNode, addLink, removeLink } from "./default_rete";
import "./widget.css";

export function useRete<T extends { destroy(): void }>(
  create: (el: HTMLElement, settings: any, data: any) => Promise<T>,
  settings: any, workgraphData: any
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
      create(container, settings, workgraphData).then((value) => {
        editorRef.current = value;
        setEditor(value);
        // Expose the editor instance to the window object for debugging
        window.editor = value;
      });
    }
  }, [container, create, workgraphData]);

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
  const [value, setValue] = useModelState<any>("value");
  const [style, setStyle] = useModelState<any>("style");
  const [settings, setSettings] = useModelState<any>("settings");
  const [ref, editor] = useRete(createEditor, settings, value);
  const model = useModel();

  React.useEffect(() => {
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
          addNode(editor.editor, msg.data)
          editor.layout(true);
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
    model.on("msg:custom", handle_custom_msg);
    return () => model.off("msg:custom", handle_custom_msg);
  }, [model, editor]);

  return (
    <div className="App">
      <div ref={ref} className="rete" style={style}></div>
    </div>
  );
});

export default { render };
