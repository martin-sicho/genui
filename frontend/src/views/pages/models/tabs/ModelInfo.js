import React from "react";
import {TaskProgressBar, TaskBadgeGroup} from "../../../../genui"
import { Col, Row } from 'reactstrap';

class ModelInfo extends React.Component {
  render() {
    const model = this.props.model;
    const tasks = this.props.tasks;

    return (<Row>
      <Col sm="12">
        <h4>Description</h4>
        <p>{model.description}</p>

        <h4>
          Tasks <TaskBadgeGroup tasks={tasks}/>
        </h4>
        <TaskProgressBar
          progressURL={this.props.apiUrls.celeryProgress}
          tasks={tasks.running}
        />
      </Col>
    </Row>)
  }
}

export default ModelInfo;