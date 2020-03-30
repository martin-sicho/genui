import React from 'react';
import { BrowserRouter, Route, Switch } from 'react-router-dom';
import DashboardLayout from './layouts/DashboardLayout';
import './vibe/scss/styles.scss';
import '../node_modules/react-grid-layout/css/styles.css'
import '../node_modules/react-resizable/css/styles.css'
import LoginPage from './views/pages/login/LoginPage';

let BACKEND_URL = new URL('http://localhost:8000');
if (process.env.REACT_APP_GENUI_BACKEND_ROOT_URL) {
  BACKEND_URL = (process.env.REACT_APP_GENUI_BACKEND_ROOT_URL);
  console.log('Using REACT_APP_GENUI_BACKEND_ROOT_URL for backend url...');
}
console.log(`Set backend URL to: ${BACKEND_URL}`);
const REMOTE_API_ROOT = new URL('api/', BACKEND_URL);
console.log(`Remote API root at: ${REMOTE_API_ROOT}`);

const generatorsURL = new URL('generators/', REMOTE_API_ROOT);
const URL_ROOTS = {
  accountsRoot: new URL('accounts/', REMOTE_API_ROOT),
  projectList : new URL('projects/', REMOTE_API_ROOT),
  compoundsRoot: new URL('compounds/', REMOTE_API_ROOT),
  compoundSetsRoot : new URL('compounds/sets/', REMOTE_API_ROOT),
  activitySetsRoot : new URL('compounds/activity/sets/', REMOTE_API_ROOT),
  qsarRoot : new URL('qsar/', REMOTE_API_ROOT),
  generatorsRoot : generatorsURL,
  drugexRoot : new URL('drugex/', generatorsURL),
  mapsRoot: new URL('maps/', REMOTE_API_ROOT),
  celeryProgress : new URL('celery-progress/', REMOTE_API_ROOT),
};

const fetchUserInfo = (callback) => {
  fetch(new URL('user/', URL_ROOTS.accountsRoot), {
    credentials: "include",
    "headers": {
      "Accept": "application/json",
    },
    "method": "GET"
  }).then(response => response.json())
    .then(data => {
      if (data.username) {
        callback(data);
      } else {
        callback(null)
      }
    })
    .catch(e => console.log(e))
};

export default function App() {
  const [user, setUser] = React.useState(null);
  const loginPagePath = '/login';
  const appPath = '/';

  return (
    <BrowserRouter>
      <Switch>
        <Route
          exact
          path={loginPagePath}
          render={
            (props) => (
              <LoginPage
                {...props}
                apiUrls={URL_ROOTS}
                fetchUserInfo={fetchUserInfo}
                setUser={setUser}
                user={user}
                loginPagePath={loginPagePath}
                appPath={appPath}
              />
            )
          }/>
        <Route path={appPath} render={
          (props) => (
            <DashboardLayout
              {...props}
              apiUrls={URL_ROOTS}
              user={user}
              setUser={setUser}
              fetchUserInfo={fetchUserInfo}
              loginPagePath={loginPagePath}
              appPath={appPath}
            />
          )
        } />
      </Switch>
    </BrowserRouter>
  );
}
