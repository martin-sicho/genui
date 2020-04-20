import React, { Component } from 'react';
import { Card, DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';
import { ResponsiveGrid } from '../../../genui/';
import { CreateNewCard } from './CreateNewCard';
import { ProjectCard } from './ProjectCard';

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

    fetchUpdates = () => {
      fetch(this.props.apiUrls.projectList, {
        credentials: "include",
        "headers": {
          "Accept": "application/json",
        },
        "method": "GET"
      })
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
                },
                credentials: "include"
            }
        ).then(response => response.json()).then(
            data => {
              let new_project = Object.assign({url : `/projects/${data.id}`}, data);
                this.setState({
                    creating: false
                });
                this.props.onProjectOpen(new_project);
                this.props.history.push(new_project.url);
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
              h : {"md" :4, "sm" : 4},
              w : {"md" : 1, "sm" : 1},
              minH : {"md" : 4, "sm" : 4},
              data : {}
      };
      // console.log(project_cards.concat(new_project_card));
      return (
      this.state.creating ? <div>Loading...</div>: <ResponsiveGrid
          items={project_cards.concat(new_project_card)}
          rowHeight={75}
          mdCols={2}
          smCols={1}
          gridID="projects-grid-layout"
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
