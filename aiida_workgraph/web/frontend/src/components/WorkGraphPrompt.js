// PromptModal.js
import 'bootstrap/dist/css/bootstrap.min.css';
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';
import { useState } from 'react';


function WorkGraphDeleteNodePrompt(props) {

  const handleConfirmClick = () => {
    props.confirmAction();
    props.setShow(false);
  };

  const handleCancelClick = () => {
    props.cancelAction();
    props.setShow(false);
  };

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
          Are you sure you want to delete node? <b>A deletion is irreversible.</b>
        </p>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={handleCancelClick}>Cancel</Button>
        <Button onClick={handleConfirmClick}>Confirm</Button>
      </Modal.Footer>
      </Modal>
  );
}
export default WorkGraphDeleteNodePrompt;
