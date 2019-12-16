import React, {Component} from 'react';
import {
    Button,
    Card,
    CardBody,
    CardFooter,
    CardHeader,
    CardSubtitle,
    Form,
    FormGroup,
    FormText,
    Input,
    Label
} from 'reactstrap';
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
              <Button color="success" onClick={() => {this.props.onProjectOpen(this.project);this.props.history.push(`/projects/${this.project.id}`)}}>Open</Button> <Button color="danger" onClick={() => this.props.onProjectDelete(this.project)}>Delete</Button>
            </CardFooter>
          </React.Fragment>
    )}
}

class Projects extends Component {
    constructor(props) {
        super(props);
        this.state = {projects : []}
    }

    componentDidMount() {
        this.fetchUpdates();
    }

    fetchUpdates = () => {
        fetch(this.props.apiUrls.projectList)
            .then(response => response.json())
            .then(this.updateProjectRoutes)
    };

    updateProjectRoutes = (data) => {
        const projects = [];
        data.forEach(
            (project) => {
              const url = '/projects/' + project.id;
              projects.push(Object.assign({url : url}, project))
            }
        );

        // this.activateProject(projects[0]);

        this.setState({
          projects : projects
        })
    };

  render() {
    return (
      <ResponsiveGrid
          items={this.state.projects}
          rowHeight={75}
      >
          {
              this.state.projects.map(item =>
                  <Card className="scrollable" key={item.id.toString()}>
                      <ProjectCard {...this.props} project={item}/>
                  </Card>
              ).concat([
                  (
                      <Card key="new-project" className="scrollable">
                        <CardHeader>Create New Project</CardHeader>
                          <CardBody>
                          <Form>
                        <FormGroup>
                          <Label for="name">Name</Label>
                          <Input type="text" name="name" id="name" placeholder="New Project Name" />
                        </FormGroup>
                        <FormGroup>
                          <Label for="description">Text Area</Label>
                          <Input type="textarea" name="description" id="description" placeholder="Description..." />
                        </FormGroup>
                        <Button>Submit</Button>
                      </Form>
                        </CardBody>
                      </Card>
                  )
              ])
          }
      </ResponsiveGrid>
    )
  }
}

export default Projects;
