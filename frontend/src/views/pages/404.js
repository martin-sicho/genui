import React from 'react';
import { NavLink } from 'react-router-dom';

const ErrorPage = () => {
  return (
    <div>
      <div className="m-t-xxl text-center">
        <h1 className="error-number">404</h1>
        <h3 className="m-b">Sorry, but we couldn't find this page...</h3>
        <NavLink to={'/'}>Go back to your projects.</NavLink>
      </div>
    </div>
  );
};

export default ErrorPage;
