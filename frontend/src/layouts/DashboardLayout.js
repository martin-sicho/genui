import React, { Component } from 'react';
import { Switch, Route } from 'react-router-dom';
import { UncontrolledDropdown, DropdownToggle, DropdownMenu, DropdownItem } from 'reactstrap';
import { Header, SidebarNav, Footer, PageContent, Avatar, PageAlert, Page, RoutedPage } from '../vibe';
import Logo from '../assets/images/vibe-logo.svg';
// import avatar1 from '../assets/images/avatar1.png';
import defaultNav from '../_nav';
import defaultRoutes from '../views';
import ContextProviders from '../vibe/components/utilities/ContextProviders';
import handleKeyAccessibility, { handleClickAccessibility } from '../vibe/helpers/handleTabAccessibility';

const MOBILE_SIZE = 992;

// TODO: this should be set during build (production vs. development)
const BACKEND_URL = new URL('http://localhost:8000/');
const REMOTE_API_ROOT = new URL('api/', BACKEND_URL);

class DashboardLayout extends Component {
  constructor(props) {
    super(props);
    this.apiUrls = {
        projectList : new URL('projects/', REMOTE_API_ROOT)
    };
    this.routes = defaultRoutes;
    this.state = {
      sidebarCollapsed: false,
      isMobile: window.innerWidth <= MOBILE_SIZE,
      showChat1: false,
      currentProject: null,
      nav: defaultNav,
    };
  }

  handleResize = () => {
    if (window.innerWidth <= MOBILE_SIZE) {
      this.setState({ sidebarCollapsed: false, isMobile: true });
    } else {
      this.setState({ isMobile: false });
    }
  };

  componentDidUpdate(prev) {
    if (this.state.isMobile && prev.location.pathname !== this.props.location.pathname) {
      this.toggleSideCollapse();
    }
  }

  componentDidMount() {
    window.addEventListener('resize', this.handleResize);
    document.addEventListener('keydown', handleKeyAccessibility);
    document.addEventListener('click', handleClickAccessibility);
  }

  activateProject = (project) => {
    const current_project = project;
    const nav = JSON.parse(JSON.stringify(defaultNav));
    nav.top.push(
      {
        name: current_project.name,
        url: current_project.url,
        icon: 'Layers',
      }
    );
    nav.top.push(
        {
          divider: true,
        },
    );
    nav.top.push(
      {
        name: "Compounds",
        url: current_project.url + "/compounds",
        icon: 'Box',
      }
    );
    nav.top.push(
      {
        name: "QSAR Models",
        url: current_project.url + "/qsar",
        icon: 'Activity',
      }
    );
    nav.top.push(
      {
        name: "Generators",
        url: current_project.url + "/generators",
        icon: 'Compass',
      }
    );
    nav.top.push(
      {
        name: "Maps",
        url: current_project.url + "/maps",
        icon: 'Map',
      }
    );

    this.setState(() => ({
        currentProject : current_project,
        nav : nav,
    }));
  };

  deleteProject = (project) => {
      return fetch(this.apiUrls.projectList + project.id + '/', {method: 'DELETE'}).then((response) => {
        if (response.ok) {
            if (this.state.currentProject && (project.id === this.state.currentProject.id)) {
              const nav = JSON.parse(JSON.stringify(defaultNav));
              this.setState(() => ({
                  currentProject : null,
                  nav : nav,
              }));
            }
        } else {
            console.log("Failed to delete project: " + project.id);
        }
      });
  };

  componentWillUnmount() {
    window.removeEventListener('resize', this.handleResize);
  }

  toggleSideCollapse = () => {
    this.setState(prevState => ({ sidebarCollapsed: !prevState.sidebarCollapsed }));
  };

  closeChat = () => {
    this.setState({ showChat1: false });
  };

  render() {
    const { sidebarCollapsed } = this.state;
    const {nav} = this.state;
    const routes = this.routes;
    const sidebarCollapsedClass = sidebarCollapsed ? 'side-menu-collapsed' : '';
    return (
      <ContextProviders>
        <div className={`app ${sidebarCollapsedClass}`}>
          <PageAlert />
          <div className="app-body">
            <SidebarNav
              nav={nav}
              logo={Logo}
              logoText="GenUI"
              isSidebarCollapsed={sidebarCollapsed}
              toggleSidebar={this.toggleSideCollapse}
              {...this.props}
            />
            <Page>
              <Header
                toggleSidebar={this.toggleSideCollapse}
                isSidebarCollapsed={sidebarCollapsed}
                routes={routes}
                {...this.props}
              >
                <HeaderNav />
              </Header>
              <PageContent>
                <Switch>
                  {routes.map(page => (
                    <Route
                        exact path={page.path}
                        key={page.key}
                        render={props => (
                            <RoutedPage
                                {...props}
                                apiUrls={this.apiUrls}
                                component={page.component}
                                title={page.name}
                                currentProject={this.state.currentProject}
                                onProjectOpen={project => this.activateProject(project)}
                                onProjectDelete={project => this.deleteProject(project)}
                            />
                        )}
                    />
                  ))}
                </Switch>
              </PageContent>
            </Page>
          </div>
          <Footer>
            <span>Copyright © 2019 Nice Dash. All rights reserved.</span>
            {/*<span>*/}
            {/*  <a href="#!">Terms</a> | <a href="#!">Privacy Policy</a>*/}
            {/*</span>*/}
            {/*<span className="ml-auto hidden-xs">*/}
            {/*  Made with{' '}*/}
            {/*  <span role="img" aria-label="taco">*/}
            {/*    🌮*/}
            {/*  </span>*/}
            {/*</span>*/}
          </Footer>
          {/*<Chat.Container>*/}
          {/*  {this.state.showChat1 && (*/}
          {/*    <Chat.ChatBox name="Messages" status="online" image={avatar1} close={this.closeChat} />*/}
          {/*  )}*/}
          {/*</Chat.Container>*/}
        </div>
      </ContextProviders>
    );
  }
}

function HeaderNav() {
  return (
    <React.Fragment>
      {/*<NavItem>*/}
      {/*  <form className="form-inline">*/}
      {/*    <input className="form-control mr-sm-1" type="search" placeholder="Search" aria-label="Search" />*/}
      {/*    <Button type="submit" className="d-none d-sm-block">*/}
      {/*      <i className="fa fa-search" />*/}
      {/*    </Button>*/}
      {/*  </form>*/}
      {/*</NavItem>*/}
      <UncontrolledDropdown nav inNavbar>
        <DropdownToggle nav caret>
          Projects
        </DropdownToggle>
        <DropdownMenu right>
          <DropdownItem>New Project</DropdownItem>
          {/*<DropdownItem>Edit</DropdownItem>*/}
          {/*<DropdownItem divider />*/}
          {/*<DropdownItem>Close</DropdownItem>*/}
          {/*<DropdownItem>Delete</DropdownItem>*/}
        </DropdownMenu>
      </UncontrolledDropdown>
      <UncontrolledDropdown nav inNavbar>
        <DropdownToggle nav>
          <Avatar size="small" color="blue" initials="JS" />
        </DropdownToggle>
        <DropdownMenu right>
          <DropdownItem>Settings</DropdownItem>
          <DropdownItem>Profile</DropdownItem>
          <DropdownItem divider />
          <DropdownItem>Log Out</DropdownItem>
        </DropdownMenu>
      </UncontrolledDropdown>
    </React.Fragment>
  );
}

export default DashboardLayout;
