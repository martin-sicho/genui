import React from 'react';
import { Badge, Col, Progress, Row } from 'reactstrap';

class MolsStats extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      molCount : 0
    }
  }

  componentDidMount() {
    fetch(this.props.moleculesURL)
      .then(response => response.json())
      .then(this.updateMolStats)
  }

  updateMolStats = (data) => {
    this.setState({molCount : data.count});
  };

  render() {
    return (
      <React.Fragment>
        <h4>Compounds</h4>
        <p>Total: {this.state.molCount}</p>
        <h4>Associated Targets</h4>
        <ul>
          {
            this.props.molset.targets.map(
              target => <li key={target.targetID}>{target.targetID}</li>
            )
          }
        </ul>
      </React.Fragment>
    )
  }
}

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
    this.interval = setInterval(this.updateProgress, 5000);
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

class MolSetTasksStatus extends React.Component {
  intervalID;

  constructor(props) {
    super(props);

    this.state = {
      tasks : null
    }
  }

  componentWillUnmount() {
    clearTimeout(this.intervalID);
  }

  componentDidMount() {
    this.updateTasks();
  }

  updateTasks = () => {
    fetch(this.props.tasksURL)
      .then(response => response.json())
      .then(data => {
        this.setState({tasks : {...data}});
        this.intervalID = setTimeout(this.updateTasks.bind(this), 5000);
      })
  };

  render() {
    const tasks = this.state.tasks;
    if (!tasks) {
      return null
    }

    const completed = [];
    const running = [];
    const errors = [];
    Object.keys(tasks).forEach(task_name => {
      tasks[task_name].forEach(task => {
        task.task_name = task_name;
        if (task.status === 'SUCCESS') {
          completed.push(task)
        } else if (['STARTED', 'RECEIVED', 'PENDING', 'RETRY', 'PROGRESS'].includes(task.status)) {
          running.push(task)
        } else if (['FAILURE', 'REVOKED'].includes(task.status)) {
          errors.push(task)
        }
      });
    });

    return (
      <React.Fragment>
        <h4>
          Tasks <TaskHeadingBadge href="#" color="primary" tasks={running}>Running</TaskHeadingBadge> <TaskHeadingBadge href="#" color="success" tasks={completed}>Completed</TaskHeadingBadge> <TaskHeadingBadge href="#" color="danger" tasks={errors}>Failed</TaskHeadingBadge>
        </h4>
        <TasksProgressOverview {...this.props} tasks={running}/>
      </React.Fragment>
    )
  }
}

class ChEMBLInfo extends React.Component {

  constructor(props) {
    super(props);

    this.molset = this.props.molset;
  }

  render() {
    return (
      <Row>
        <Col sm="12">
          <h4>Description</h4>
          <p>{this.molset.description}</p>
          <MolsStats molset={this.props.molset} moleculesURL={this.props.moleculesURL} />
          <MolSetTasksStatus progressURL={this.props.apiUrls.celeryProgress} molset={this.props.molset} tasksURL={this.props.tasksURL}/>
        </Col>
      </Row>
    );
  }
}

export default ChEMBLInfo;