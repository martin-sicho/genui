// import Dashboard from './pages/Dashboard';
// import Buttons from './elements/Buttons';
// import Alerts from './elements/Alerts';
// import Grid from './elements/Grid';
// import Typography from './elements/Typography';
// import Cards from './elements/Cards';
// import Tabs from './elements/Tabs';
// import Tables from './elements/Tables';
// import Breadcrumbs from './elements/Breadcrumbs';
// import Forms from './elements/Forms';
// import Loaders from './elements/Loaders';
// import Avatars from './elements/Avatars';
// import Invoice from './pages/Invoice';
// import Analytics from './pages/Analytics';
// import CmsPage from './pages/Cms';
// import Widgets from './pages/Widgets';
// import BlankPage from './pages/BlankPage';
// import SubNav from './pages/SubNav';
// import Feed from './pages/Feed';
// import Modals from './elements/Modals';
// import ProgressBars from './elements/ProgressBars';
// import PaginationPage from './elements/Pagination';
import Projects from "./pages/projects/Projects";
import ProjectOverview from "./pages/projects/ProjectOverview";
import ErrorPage from "./pages/404";
import Compounds from './pages/compounds/Compounds';
import Models from "./pages/models/Models"
import DrugExPage from './pages/generators/drugex/DrugEx';
import MapCreator from './pages/maps/create/MapCreator';
import Maps from './pages/maps/display/MapDashboard';

// See React Router documentation for details: https://reacttraining.com/react-router/web/api/Route
const pageList = [
  {
    name: 'Projects',
    path: ["/", "/projects"],
    key: 'projects',
    component: Projects,
  },
  {
    name: 'Project Overview',
    path: ["/projects/:project"],
    key: 'projects-name',
    component: ProjectOverview,
  },
  {
    name: 'Compounds',
    path: ["/projects/:project/compounds"],
    key: 'data',
    component: Compounds,
  },
  {
    name: 'QSAR Models',
    path: ["/projects/:project/qsar"],
    key: 'qsar',
    component: Models,
  },
  {
    name: 'DrugEx',
    path: ["/projects/:project/generators/drugex"],
    key: 'generators',
    component: DrugExPage,
  },
  {
    name: 'Map Creator',
    path: ["/projects/:project/maps/creator"],
    key: 'maps-creator',
    component: MapCreator,
  },
  {
    name: 'Map Explorer',
    path: ["/projects/:project/maps/explorer"],
    key: 'maps-explorer',
    component: Maps,
  },
  {
    name: '404',
    key: 'NotFound-404',
    component: ErrorPage,
  },
  // {
  //   name: 'Buttons',
  //   path: '/elements/buttons',
  //   component: Buttons,
  // },
  // {
  //   name: 'Alerts',
  //   path: '/elements/alerts',
  //   component: Alerts,
  // },
  // {
  //   name: 'Grid',
  //   path: '/elements/grid',
  //   component: Grid,
  // },
  // {
  //   name: 'Typography',
  //   path: '/elements/typography',
  //   component: Typography,
  // },
  // {
  //   name: 'Cards',
  //   path: '/elements/cards',
  //   component: Cards,
  // },
  // {
  //   name: 'Tabs',
  //   path: '/elements/summaries',
  //   component: Tabs,
  // },
  // {
  //   name: 'Tables',
  //   path: '/elements/tables',
  //   component: Tables,
  // },
  // {
  //   name: 'Progress Bars',
  //   path: '/elements/progressbars',
  //   component: ProgressBars,
  // },
  // {
  //   name: 'Pagination',
  //   path: '/elements/pagination',
  //   component: PaginationPage,
  // },
  // {
  //   name: 'Modals',
  //   path: '/elements/modals',
  //   component: Modals,
  // },
  // {
  //   name: 'Breadcrumbs',
  //   path: '/elements/breadcrumbs',
  //   component: Breadcrumbs,
  // },
  // {
  //   name: 'Forms',
  //   path: '/elements/forms',
  //   component: Forms,
  // },
  // {
  //   name: 'Loaders',
  //   path: '/elements/loaders',
  //   component: Loaders,
  // },
  // {
  //   name: 'Avatars',
  //   path: '/elements/avatars',
  //   component: Avatars,
  // },
  // {
  //   name: 'Blank',
  //   path: '/pages/blank',
  //   component: BlankPage,
  // },
  // {
  //   name: 'Sub Navigation',
  //   path: '/pages/subnav',
  //   component: SubNav,
  // },
  // {
  //   name: 'Analytics',
  //   path: '/apps/analytics',
  //   component: Analytics,
  // },
  // {
  //   name: 'Activity Feed',
  //   path: '/apps/feed',
  //   component: Feed,
  // },
  // {
  //   name: 'Invoice',
  //   path: '/apps/invoice',
  //   component: Invoice,
  // },
  // {
  //   name: 'CMS',
  //   path: '/apps/cms',
  //   component: CmsPage,
  // },
  // {
  //   name: 'Widgets',
  //   path: '/widgets',
  //   component: Widgets,
  // },
];

export default pageList;
