import React, { Component } from 'react';
import {
  Row,
  Button,
  Col,
  Card,
  CardHeader,
  CardFooter,
  CardBody,
  CardSubtitle
} from 'reactstrap';

const KEY = require('weak-key');

function ProjectCard(props) {
    const project = props.project;
    const created = new Date(project.created);
    const updated = new Date(project.updated);

    return (
        <Card>
            <CardHeader>{project.name}</CardHeader>
            <CardBody>
                <CardSubtitle>
                    <p>
                        Created: {
                        created.toLocaleDateString()
                        + ' – ' + created.toLocaleTimeString()
                    }
                        <br/>
                        Last Update: {
                        updated.toLocaleDateString()
                        + ' – ' + updated.toLocaleTimeString()
                    }
                    </p>
                </CardSubtitle>
                <p>
                    {project.description}
                </p>
            </CardBody>
            <CardFooter>
              <Button color="success">Open</Button> <Button color="primary">Edit</Button> <Button color="secondary">Delete</Button>
            </CardFooter>
          </Card>
    );
}

function ProjectCol(props) {
    return (
        <Col md={props.colWidth}>
            {
                props.projects.map(project =>
                    <ProjectCard key={project.id} project={project}/>
                )
            }
        </Col>
    );
}

function ProjectGrid(props) {
    const col_width = 12 / props.cardsPerCol;
    const projects = props.projects;
    // const active_project = props.currentProject; // TODO: highlight this somehow
    let current_batch = [];
    let cols = [];
    projects.map((project, index) => {
        current_batch.push(project);
        if (index+1 % props.cardsPerCol) {
            cols.push(<ProjectCol colWidth={col_width} projects={current_batch}/>);
            current_batch = [];
        }
    });

    return (
        cols.map(col => <React.Fragment key={KEY(col)}>{col}</React.Fragment>)
    );
}

class Projects extends Component {
  constructor(props) {
    super(props);
    this.cardsPerCol = 2;
  }

  render() {
    return (
      <Row>
          <ProjectGrid
              cardsPerCol={this.cardsPerCol}
              projects={this.props.projects}
              currentProject={this.props.currentProject}
          />
      </Row>
    );
  }
}

export default Projects;
