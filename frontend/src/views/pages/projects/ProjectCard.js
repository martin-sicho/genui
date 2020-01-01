import React from 'react';
import { Button, CardBody, CardFooter, CardHeader, CardSubtitle } from 'reactstrap';

export class ProjectCard extends React.Component {

  constructor(props) {
    super(props);
    this.project = props.project;
    this.created = new Date(this.project.created);
    this.updated = new Date(this.project.updated);
  }

  render() {
    return (
      <React.Fragment>
        <CardHeader>{this.project.name}</CardHeader>
        <CardBody className="scrollable">
          <CardSubtitle>
            <p>
              Created: {
              this.created.toLocaleDateString()
              + ' – ' + this.created.toLocaleTimeString()
            }
              <br/>
              Last Update: {
              this.updated.toLocaleDateString()
              + ' – ' + this.updated.toLocaleTimeString()
            }
            </p>
          </CardSubtitle>
          <p>
            {this.project.description}
          </p>
        </CardBody>
        <CardFooter>
          <Button color="success" onClick={() => {
            this.props.onProjectOpen(this.project);
            this.props.history.push(`/projects/${this.project.id}`);
          }}>Open</Button> <Button color="danger"
                                   onClick={() => this.props.onProjectDelete(this.project)}>Delete</Button>
        </CardFooter>
      </React.Fragment>
    );
  }
}