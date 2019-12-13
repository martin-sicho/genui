import React, {Component} from 'react';
import {Button, Card, CardBody, CardFooter, CardHeader, CardSubtitle} from 'reactstrap';
import {ResponsiveGrid} from "../../vibe/components/grid/ResponsiveGrid";
import "./Projects.css";

class ProjectCard extends React.Component {

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
            <CardBody>
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
              <Button color="success">Open</Button> <Button color="primary">Edit</Button> <Button color="secondary">Delete</Button>
            </CardFooter>
          </React.Fragment>
    )}
}

class Projects extends Component {

  render() {
    return (
      <ResponsiveGrid
          items={this.props.projects}
          rowHeight={75}
      >
          {
              this.props.projects.map(item =>
                  <Card className="scrollable" key={item.id.toString()}>
                      <ProjectCard project={item}/>
                  </Card>
              )
          }
      </ResponsiveGrid>
    )
  }
}

export default Projects;
