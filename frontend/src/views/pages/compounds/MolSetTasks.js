import React from 'react';
import { Badge, Progress } from 'reactstrap';

function TaskHeadingBadge(props) {
  if (props.tasks.length === 0) return null;

  return (
    <Badge {...props}>{props.tasks.length} {props.children}</Badge>
  )
}

class TasksProgressOverview extends React.Component {
  constructor(props) {
    super(props);

    this.state = {
      progressData : []
    }
  }

  componentDidMount() {
    this.updateProgress();
    this.interval = setInterval(this.updateProgress, 2000);
  }

  componentWillUnmount() {
    clearInterval(this.interval);
  }

  updateProgress = () => {
    const tasks = this.props.tasks;
    const progressData = [];
    tasks.forEach(task => {
      const url = new URL(task.task_id + '/', this.props.progressURL);
      fetch(url)
        .then(response => response.json())
        .then(data => {
          data.task = task;
          progressData.push(data);

          // FIXME: this set state should not happen if the component is unmounted -> the fetch needs to be cancelled properly
          this.setState(state => {
            return {
              progressData : progressData
            };
          });
        })
      ;
    });
  };

  render() {
    const progress = this.state.progressData;
    progress.sort((a, b) => (a.task.task_id > b.task.task_id) ? 1 : -1);

    return (
      <React.Fragment>
        {
          progress.map(data => (
            <React.Fragment key={data.task.task_id}>
              <div className="text-center">{data.task.task_name} ({data.progress.percent}%)</div>
              <Progress value={data.progress.percent} />
            </React.Fragment>
          ))
        }
      </React.Fragment>
    )
  }
}

class MolSetTasks extends React.Component {

  render() {
    const tasks = this.props.tasks;
    if (!tasks) {
      return null
    }

    return (
      <React.Fragment>
        <h4>
          Tasks <TaskHeadingBadge href="#" color="primary" tasks={tasks.running}>Running</TaskHeadingBadge> <TaskHeadingBadge href="#" color="success" tasks={tasks.completed}>Completed</TaskHeadingBadge> <TaskHeadingBadge href="#" color="danger" tasks={tasks.errors}>Failed</TaskHeadingBadge>
        </h4>
        {
          tasks.running.length > 0 ?
            <TasksProgressOverview {...this.props} tasks={tasks.running}/>
            : null
        }
      </React.Fragment>
    )
  }
}

export default MolSetTasks;