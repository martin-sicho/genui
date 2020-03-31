import React, { Component } from 'react';
import { Switch, Route, Redirect, useHistory } from 'react-router-dom';
import { UncontrolledDropdown, DropdownToggle, DropdownMenu, DropdownItem } from 'reactstrap';
import { Header, SidebarNav, Footer, PageContent, Avatar, PageAlert, Page} from '../vibe';
import {RoutedPage} from '../genui/'
import Logo from '../assets/images/vibe-logo.svg';
// import avatar1 from '../assets/images/avatar1.png';
import defaultNav from '../_nav';
import defaultRoutes from '../views';
import ContextProviders from '../vibe/components/utilities/ContextProviders';
import handleKeyAccessibility, { handleClickAccessibility } from '../vibe/helpers/handleTabAccessibility';
import LogInManager from '../views/pages/login/LoginManager';

const MOBILE_SIZE = 992;

class DashboardLayout extends Component {
  constructor(props) {
    super(props);
    this.apiUrls = this.props.apiUrls;
    this.routes = defaultRoutes;
    this.state = {
      pageTitle : "GenUI",
      sidebarCollapsed: false,
      isMobile: window.innerWidth <= MOBILE_SIZE,
      showChat1: false,
      currentProject: null,
      nav: defaultNav,
      headerComponent: null,
      failedToLogIn: false,
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
      this.handleResize();
    }
  }

  componentDidMount() {
    // window.addEventListener('resize', this.handleResize);
    document.addEventListener('keydown', handleKeyAccessibility);
    document.addEventListener('click', handleClickAccessibility);
    if (!this.props.user) {
      this.props.fetchUserInfo((userData) => {
        if (userData) {
          this.props.setUser(userData);
        } else {
          this.setState({failedToLogIn: true})
        }
      })
    }
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
        url: current_project.url + "compounds",
        icon: 'Box',
      }
    );
    nav.top.push(
      {
        name: "QSAR Models",
        url: current_project.url + "qsar",
        icon: 'Activity',
      }
    );
    nav.top.push(
      {
        name: "Generators",
        icon: 'Compass',
        children: [
          {
            name: 'DrugEx',
            url: current_project.url + "generators/drugex",
          }
        ],
      },
    );
    nav.top.push(
      {
        name: "Maps",
        icon: 'Map',
        children: [
          {
            name: 'Creator',
            url: current_project.url + "maps/creator",
          },
          {
            name: 'Explorer',
            url: current_project.url + "maps/explorer",
          }
        ],
      }
    );

    this.setState(() => ({
        currentProject : current_project,
        nav : nav,
    }));
  };

  deleteProject = (project) => {
      return fetch(this.apiUrls.projectList + project.id + '/', {method: 'DELETE', credentials: "include"}).then((response) => {
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
    if (this.state.failedToLogIn) {
      return <Redirect to={this.props.loginPagePath}/>
    }

    if (!this.props.user) {
      return <div/>
    }

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
                <HeaderNav {...this.props} injected={this.injectContentToHeader} />
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
                                handlePageTitleChange={this.handlePageTitleChange}
                                apiUrls={this.apiUrls}
                                component={page.component}
                                title={page.name}
                                currentProject={this.state.currentProject}
                                onProjectOpen={project => this.activateProject(project)}
                                onProjectDelete={project => this.deleteProject(project)}
                                onHeaderChange={this.handleHeaderChange}
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

    injectContentToHeader = () => {
        return this.state.headerComponent
    };

    handleHeaderChange = (component) => {
      this.setState({
          headerComponent : component
      })
    };

    handlePageTitleChange = (newTitle) => {
      document.title = `GenUI > ${newTitle}`;
      this.setState(() => ({pageTitle : newTitle}));
    }
}

function HeaderNav(props) {
   const history = useHistory();
   const Injected = props.injected ? props.injected : React.Fragment;
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
      <Injected/>
      <UncontrolledDropdown nav inNavbar>
        <DropdownToggle nav>
          <Avatar size="small" color="blue" initials="JS" />
        </DropdownToggle>
        <DropdownMenu right>
          {/*<DropdownItem>Settings</DropdownItem>*/}
          <DropdownItem onClick={(e) => {
            e.preventDefault();
            history.push(props.loginPagePath);
          }}>Profile</DropdownItem>
          <DropdownItem divider />
          <LogInManager
            {...props}
          >
            {
              (sendLogInRequest, sendLogOutRequest) => (
                <DropdownItem onClick={(e) => {
                  e.preventDefault();
                  sendLogOutRequest();
                  history.push(props.loginPagePath)
                }}>Log Out</DropdownItem>
              )
            }
          </LogInManager>
        </DropdownMenu>
      </UncontrolledDropdown>
    </React.Fragment>
  );
}

export default DashboardLayout;
