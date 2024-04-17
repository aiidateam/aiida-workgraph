import * as React from "react";
import { Button, Switch } from "antd";
import { createRender, useModel, useModelState } from "@anywidget/react";
import { createEditor, loadJSON, removeNode, addNode, addLink, removeLink } from "./default_rete";
import "./widget.css";

export function useRete<T extends { destroy(): void }>(
  create: (el: HTMLElement) => Promise<T>,
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
      create(container).then((value) => {
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

const render = createRender(() => {
  const [value, setValue] = useModelState<any>("value");
  const [ref, editor, isInitialized] = useRete(createEditor);
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
          const node = editor.editor.getNode(context.data.id)
          value.nodes[node.label].position = [view.position.x, view.position.y]
          value.uuid = Math.random().toString(36).substring(7)
          editor.editor.uuid = value.uuid
          console.log("value: ", value)
          // copy the value using json parse and stringify to trigger a re-render
          setValue(JSON.parse(JSON.stringify(value)))
        }
      }
      return context
    })
  }, [isInitialized]);

  React.useEffect(() => {
    model.on("msg:custom", handle_custom_msg);
    return () => model.off("msg:custom", handle_custom_msg);
  }, [model, editor]);

  React.useEffect(() => {
    console.log("Value changed", value);
    if (!editor) {
      console.log("Editor not ready");
      return;
    };
    // same uuid skip
    if (editor.editor.uuid === value.uuid) {
      console.log("Same uuid, skip");
      return;
    }
    loadJSON(editor.editor, editor.area, editor.layout, value);
  }, [isInitialized, value]);

  return (
    <div className="App">
      <div>
        <Button onClick={() => editor?.layout(true)}>Arrange</Button>
      </div>
      <div ref={ref} className="rete"></div>
    </div>
  );
});

export default { render };
