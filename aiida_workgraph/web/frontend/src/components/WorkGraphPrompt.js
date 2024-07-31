// PromptModal.js
import 'bootstrap/dist/css/bootstrap.min.css';
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';


function WorkGraphDeleteNodePrompt(props) {
  return (
      <Modal 
        {...props}
      >
      <Modal.Header closeButton>
        <Modal.Title id="contained-modal-title-vcenter">
          Confirm deletion 
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <p>
          Are you sure you want to delete node {props.item.pk}? A deletion is irreversible.
        </p>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={props.onNoClick}>No</Button>
        <Button onClick={props.onYesClick}>Yes</Button>
      </Modal.Footer>
      </Modal>
  );
}
//function WorkGraphDeleteNodePrompt({ openOnInit }) {
//    const [open, setOpen] = useState(false);
//    setOpen(openOnInit)
//
//    const handleYesClick = () => {
//      setOpen(false)
//    };
//    const handleNoClick = () => {
//      setOpen(true);
//    };
//
//    return ( 
//      <div
//        className="modal show"
//        style={{ display: 'block', position: 'initial' }}
//      >
//        <Modal isOpen={open}>
//            <h2>Are you sure you want to delete node ...</h2>
//            <button onClick={handleYesClick}>Yes</button>
//            <button onClick={handleNoClick}>No</button>
//        </Modal>
//      </div>);
//}

export default WorkGraphDeleteNodePrompt;
