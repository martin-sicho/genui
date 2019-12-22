import React, {Component} from 'react';
import {
    Button,
    Card,
    CardBody,
    CardFooter,
    CardHeader,
    CardSubtitle, DropdownItem, DropdownMenu, DropdownToggle,
    Form,
    FormGroup,
    Input,
    Label, UncontrolledAlert, UncontrolledDropdown
} from 'reactstrap';
import {ResponsiveGrid} from "../../vibe/components/grid/ResponsiveGrid";
import "./Projects.css";
import { useFormik } from 'formik';
import * as Yup from 'yup';

function HeaderNav(props) {
    return (<UncontrolledDropdown nav inNavbar>
        <DropdownToggle nav caret>
          Actions
        </DropdownToggle>
        <DropdownMenu right>
          <DropdownItem onClick={() => document.getElementById("new-proj-card").scrollIntoView()}>New Project</DropdownItem>
          <DropdownItem divider />
            <UncontrolledDropdown>
                <DropdownToggle nav>Open...</DropdownToggle>
                <DropdownMenu>
                    {
                        props.projects.map(project =>
                            (<DropdownItem
                                key={project.id}
                                onClick={() => {props.onProjectOpen(project);props.history.push(`/projects/${project.id}`)}}
                            >
                                {project.name}
                            </DropdownItem>)
                        )
                    }
                </DropdownMenu>
            </UncontrolledDropdown>
        </DropdownMenu>
      </UncontrolledDropdown>)
}

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
            , isLoading : true
        }
    }

    componentDidMount() {
        this.fetchUpdates();
    }

    componentWillUnmount() {
        this.props.onHeaderChange(null);
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
              const url = '/projects/' + project.id + '/';
              projects.push(Object.assign({url : url}, project))
            }
        );

        // this.activateProject(projects[0]);

        this.setState({
            projects : projects,
            isLoading : false
        });
        this.props.onHeaderChange(<HeaderNav {...this.props} projects={projects}/>);
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
      if (this.state.isLoading) {
          return <div>Loading...</div>
      }

      const project_cards = this.state.projects.map(project => ({
              id : project.id,
              h : {"md" : 3, "sm" : 3},
              w : {"md" : 1, "sm" : 1},
              minH : {"md" : 3, "sm" : 3},
              data : project
          }));
      const new_project_card = {
              id : "new-project",
              h : {"md" : 3, "sm" : 3},
              w : {"md" : 1, "sm" : 1},
              minH : {"md" : 3, "sm" : 3},
              data : {}
      };
      // console.log(project_cards.concat(new_project_card));
      return (
      this.state.creating ? <div>Loading...</div>: <ResponsiveGrid
          items={project_cards.concat(new_project_card)}
          rowHeight={100}
          mdCols={2}
          smCols={1}
      >
          {
              project_cards.map(item =>
                  <Card key={item.id.toString()}>
                      <ProjectCard {...this.props} project={item.data} onProjectDelete={project => {this.props.onProjectDelete(project).then(this.fetchUpdates)}}/>
                  </Card>
              ).concat([
                  (
                      <Card key="new-project" id="new-proj-card">
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
