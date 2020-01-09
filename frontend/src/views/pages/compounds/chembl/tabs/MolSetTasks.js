import React from 'react';
import {TaskBadgeGroup, TaskProgressBar} from '../../../../../genui'

class MolSetTasks extends React.Component {

  render() {
    const tasks = this.props.tasks;
    if (!tasks) {
      return null
    }

    return (
      <React.Fragment>
        <h4>
          Tasks <TaskBadgeGroup tasks={tasks}/>
        </h4>
        <TaskProgressBar
          progressURL={this.props.progressURL}
          tasks={tasks.running}
        />
      </React.Fragment>
    )
  }
}

export default MolSetTasks;