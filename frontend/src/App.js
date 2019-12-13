import React from 'react';
import { BrowserRouter, Route, Switch } from 'react-router-dom';
import DashboardLayout from './layouts/DashboardLayout';
import './vibe/scss/styles.scss';
import '../node_modules/react-grid-layout/css/styles.css'
import '../node_modules/react-resizable/css/styles.css'

export default function App() {
  return (
    <BrowserRouter>
      <Switch>
        <Route component={DashboardLayout} />
      </Switch>
    </BrowserRouter>
  );
}
