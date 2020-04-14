import {
  Badge,
  Button,
  Card,
  CardBody,
  Modal,
  ModalBody,
  ModalFooter,
  ModalHeader,
  UncontrolledCollapse,
} from 'reactstrap';
import React from 'react';

function TaskBadge(props) {
  const [open, setOpen] = React.useState(false);
  const toggle = e => {e.preventDefault();setOpen(!open)};

  if (props.tasks.length === 0) return null;

  return (
    <React.Fragment>
      <Badge {...props} onClick={toggle}>{props.tasks.length} {props.children}</Badge>
      <Modal isOpen={open} toggle={toggle}>
        <ModalHeader toggle={toggle}>Tasks</ModalHeader>
        <ModalBody>
          {
            props.tasks.map((task, index) => {
              return (
                <div key={task.task_id}>
                  <Button color="primary" id={`toggler-${index}`} style={{ marginBottom: '1rem' }}>
                    {task.task_name}
                  </Button>
                  <UncontrolledCollapse toggler={`#toggler-${index}`}>
                    <Card>
                      <CardBody>
                        <ul>
                          <li>
                            <strong>ID:</strong> {task.task_id}
                          </li>
                          <li>
                            <strong>Status:</strong> {task.status}
                          </li>
                          <li>
                            <strong>Result:</strong> {task.result}
                          </li>
                          {task.traceback ? (
                            <li>
                              <strong>Traceback:</strong> {task.traceback}
                            </li>
                          ) : null}
                        </ul>
                      </CardBody>
                    </Card>
                  </UncontrolledCollapse>
                </div>
              )
            })
          }
        </ModalBody>
        <ModalFooter>
          <Button color="primary" onClick={toggle}>Close</Button>{' '}
          {/*<Button color="secondary" onClick={toggle}>Cancel</Button>*/}
        </ModalFooter>
      </Modal>
    </React.Fragment>
  )
}

export default TaskBadge;