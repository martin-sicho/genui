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
    Input,
    Label, UncontrolledAlert
} from 'reactstrap';
import {ResponsiveGrid} from "../../vibe/components/grid/ResponsiveGrid";
import "./Projects.css";
import { useFormik } from 'formik';
import * as Yup from 'yup';

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

function CreateNewForm(props) {
    const formik = useFormik({
        initialValues: {
            name: 'New Project Name',
            description: 'Short description of the project...'
        },
        validationSchema: Yup.object({
          name: Yup.string()
              .max(256, 'Name must be less than 256 character long.')
              .required('Name is required.'),
          description: Yup.string()
              .max(10000, 'Description must be 10,000 characters or less')
        }),
        onSubmit: values => {
            props.handleCreate(values);
        },
    });

    return (
        <Form onSubmit={formik.handleSubmit} className="unDraggable">
            <FormGroup>
                <Label htmlFor="name">Name</Label>
                <Input
                    {...formik.getFieldProps('name')}
                    type="text"
                    />
            </FormGroup>
            {formik.touched.name && formik.errors.name ? <UncontrolledAlert color="danger">{formik.errors.name}</UncontrolledAlert> : null}
            <FormGroup>
                <Label htmlFor="description">Description</Label>
                <Input
                    {...formik.getFieldProps('description')}
                    type="textarea"
                    />
            </FormGroup>
            {formik.touched.description && formik.errors.description ? <UncontrolledAlert color="danger">{formik.errors.description}</UncontrolledAlert> : null}
            <Button type="submit" color="primary">Create & Open</Button>
        </Form>
    );
}

class CreateNewCard extends React.Component {

    render() {
        return (
            <React.Fragment>
                <CardHeader>Create New Project</CardHeader>
                <CardBody>
                    <CreateNewForm {...this.props}/>
                </CardBody>
            </React.Fragment>
        );
    }
}

class Projects extends Component {
    constructor(props) {
        super(props);
        this.state = {
            projects : []
            , creating : false
        }
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

    handleCreate = (values) => {
        this.setState({ creating: true });
        fetch(
            this.props.apiUrls.projectList
            , {
                method: 'POST'
                , body: JSON.stringify(values)
                , headers: {
                  'Content-Type': 'application/json'
                }
            }
        ).then(response => response.json()).then(
            data => {
                let new_project = Object.assign({url : `/projects/${data.id}`}, data);
                this.setState({
                    creating: false
                });
                this.props.onProjectOpen(new_project);
                this.props.history.push(new_project.url)
            }
        )
        ;
    };

  render() {
    return (
      this.state.creating ? <div>Loading...</div>: <ResponsiveGrid
          items={this.state.projects}
          rowHeight={75}
      >
          {
              this.state.projects.map(item =>
                  <Card className="scrollable" key={item.id.toString()}>
                      <ProjectCard {...this.props} project={item} onProjectDelete={project => {this.props.onProjectDelete(project).then(this.fetchUpdates)}}/>
                  </Card>
              ).concat([
                  (
                      <Card key="new-project" className="scrollable">
                        <CreateNewCard handleCreate={this.handleCreate} />
                      </Card>
                  )
              ])
          }
      </ResponsiveGrid>
    )
  }
}

export default Projects;
