/*
 * Class:        JDBCManager
 * Description:  Interface to access databases
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.util;
import java.io.*;
import java.net.URL;
import java.sql.*;
import javax.sql.DataSource;
import java.util.Properties;
import javax.naming.*;

/**
 * This class provides some facilities to connect to a SQL database and to
 * retrieve data stored in it. JDBC provides a standardized interface for
 * accessing a database independently of a specific database management
 * system (DBMS). The user of JDBC must create a  Connection object used to
 * send SQL queries to the underlying DBMS, but the creation of the
 * connection adds a DBMS-specific portion in the application. This class
 * helps the developer in moving the DBMS-specific information out of the
 * source code by storing it in a properties file. The methods in this class
 * can read such a properties file and establish the JDBC connection. The
 * connection can be made by using a  DataSource obtained through a JNDI
 * server, or by a JDBC URI associated with a driver class. Therefore, the
 * properties used to connect to the database must be a JNDI name
 * (<tt>jdbc.jndi-name</tt>), or a driver to load (<tt>jdbc.driver</tt>) with
 * the URI of a database (<tt>jdbc.uri</tt>).
 *
 * @code
 *
 *    jdbc.driver=com.mysql.jdbc.Driver\\
 * jdbc.uri=jdbc:mysql://mysql.iro.umontreal.ca/database?user=foo&password=bar
 *
 * @endcode
 *
 * The connection is established using the  #connectToDatabase(Properties)
 * method. Shortcut methods are also available to read the properties from a
 * file or a resource before establishing the connection. This class also
 * provides shortcut methods to read data from a database and to copy the
 * data into Java arrays.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class JDBCManager {

   /**
    * Connects to the database using the properties `prop` and returns the
    * an object representing the connection. The properties stored in
    * `prop` must be a JNDI name (<tt>jdbc.jndi-name</tt>), or the name of
    * a driver (<tt>jdbc.driver</tt>) to load and the URI of the database
    * (<tt>jdbc.uri</tt>). When a JNDI name is given, this method
    * constructs a context using the nullary constructor of
    * InitialContext, uses the context to get a  DataSource object, and
    * uses the data source to obtain a connection. This method assumes
    * that JNDI is configured correctly; see the class  InitialContext for
    * more information about configuring JNDI. If no JNDI name is
    * specified, the method looks for a JDBC URI. If a driver class name
    * is specified along with the URI, the corresponding driver is loaded
    * and registered with the JDBC  DriverManager. The driver manager is
    * then used to obtain the connection using the URI. This method throws
    * an  SQLException if the connection failed and an
    * IllegalArgumentException if the properties do not contain the
    * required values.
    *  @param prop         the properties to connect to the database.
    *  @return the connection to the database.
    *
    *  @exception SQLException if the connection failed.
    *  @exception IllegalArgumentException if the properties do not
    * contain the require values.
    */
   public static Connection connectToDatabase (Properties prop)
            throws SQLException {
      Connection connection = null;
      String jndiName;
        
      if ((jndiName = prop.getProperty ("jdbc.jndi-name")) != null)
      {
         try
         {
            InitialContext context = new InitialContext();
            connection = ((DataSource) context.lookup (jndiName)).getConnection();
         }
         catch (NamingException e)
         {
            throw new IllegalArgumentException
                  ("The jdbc.jndi-name property refers to the invalid name " + jndiName);
         }
      }
      else
      {
         String driver = prop.getProperty ("jdbc.driver");
         String uri = prop.getProperty ("jdbc.uri");
         if (uri != null)
         {
            if (driver != null) {
               try
               {
                  Class driverClass = Class.forName (driver);
                  if (!Driver.class.isAssignableFrom (driverClass))
                     throw new IllegalArgumentException
                        ("The driver name " + driver +
                     " does not correspond to a class implementing the java.sql.Driver interface");
                  // Needed by some buggy drivers
                  driverClass.newInstance();
               }
               catch (ClassNotFoundException cnfe) {
                  throw new IllegalArgumentException ("Could not find the driver class " + driver);
               }
               catch (IllegalAccessException iae) {
                  throw new IllegalArgumentException
                        ("An illegal access prevented the instantiation of driver class " + driver);
               }
               catch (InstantiationException ie) {
                  throw new IllegalArgumentException
                        ("An instantiation exception prevented the instantiation of driver class " + driver +
                                ": " + ie.getMessage());
               }
            }
            connection = DriverManager.getConnection (uri);
         }
         else
         {
            throw new IllegalArgumentException
                ("The jdbc.driver and jdbc.uri properties must be given if jdbc.jndi-name is not set");
         }            
      }

      return connection;
   }

   /**
    * Returns a connection to the database using the properties read from
    * stream `is`. This method loads the properties from the given stream,
    * and calls  #connectToDatabase(Properties) to establish the
    * connection.
    *  @param is           the stream to read for the properties.
    *  @return the connection to the database.
    *
    *  @exception SQLException if the connection failed.
    *  @exception IOException if the stream can not be read correctly.
    *  @exception IllegalArgumentException if the properties do not
    * contain the require values.
    */
   public static Connection connectToDatabase (InputStream is)
            throws IOException, SQLException {
      Properties prop = new Properties();
        
      prop.load (is);
        
      return connectToDatabase (prop);
   }

   /**
    * Equivalent to  {@link #connectToDatabase(InputStream)
    * connectToDatabase(url.openStream())}.
    */
   public static Connection connectToDatabase (URL url)
            throws IOException, SQLException {
    
      InputStream is = url.openStream();
      try {
         return connectToDatabase (is);
      }
      finally {
         is.close();
      }
   }

   /**
    * Equivalent to  {@link #connectToDatabase(InputStream)
    * connectToDatabase(new FileInputStream (file))}.
    */
   public static Connection connectToDatabase (File file)
            throws IOException, SQLException {
    
      FileInputStream is = new FileInputStream (file);
      try {
         return connectToDatabase (is);
      }
      finally {
         is.close();
      }
   }

   /**
    * Equivalent to  {@link #connectToDatabase(InputStream)
    * connectToDatabase(new FileInputStream (fileName))}.
    */
   public static Connection connectToDatabase (String fileName)
            throws IOException, SQLException {

      FileInputStream is = new FileInputStream (fileName);
      try {
         return connectToDatabase (is);
      }
      finally {
         is.close();
      }
   }

   /**
    * Uses  #connectToDatabase(InputStream) with the stream obtained from
    * the resource `resource`. This method searches the file `resource` on
    * the class path, opens the first resource found, and extracts
    * properties from it. It then uses  #connectToDatabase(Properties) to
    * establish the connection.
    */
   public static Connection connectToDatabaseFromResource (String resource)
            throws IOException, SQLException {
        
      InputStream is = JDBCManager.class.getClassLoader().getResourceAsStream (resource);
      try {
         return connectToDatabase (is);
      }
      finally {
         is.close();
      }
   }

   /**
    * Copies the result of the SQL query `query` into an array of
    * double-precision values. This method uses the statement `stmt` to
    * execute the given query, and assumes that the first column of the
    * result set contains double-precision values. Each row of the result
    * set then becomes an element of an array of double-precision values
    * which is returned by this method. This method throws an
    * SQLException if the query is not valid.
    *  @param stmt         the statement used to make the query.
    *  @param query        the SQL query to execute.
    *  @return the first column of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static double[] readDoubleData (Statement stmt, String query)
            throws SQLException {
      ResultSet rs = stmt.executeQuery (query);
      rs.last();
      double[] res = new double[rs.getRow()];
      rs.first();
        
      for (int i = 0; i < res.length; i++)
      {
         res[i] = rs.getDouble (1);
         rs.next();
      }
      rs.close();
            
      return res;
   }

   /**
    * Copies the result of the SQL query `query` into an array of
    * double-precision values. This method uses the active connection
    * `connection` to create a statement, and passes this statement, with
    * the query, to  #readDoubleData(Statement,String), which returns an
    * array of double-precision values.
    *  @param connection   the active connection to the database.
    *  @param query        the SQL query to execute.
    *  @return the first column of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static double[] readDoubleData (Connection connection,
                                          String query)
            throws SQLException {
      Statement stmt = connection.createStatement
      (ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
      try {
         return readDoubleData (stmt, query);
      }
      finally {
         stmt.close();
      }
   }

   /**
    * Returns the values of the column `column` of the table `table`. This
    * method is equivalent to  {@link #readDoubleData(Statement,String)
    * readDoubleData(stmt, "SELECT column FROM table")}.
    */
   public static double[] readDoubleData (Statement stmt, String table,
                                          String column)
            throws SQLException {
        final String query = "SELECT " + column + " FROM " + table;
        
        return readDoubleData (stmt, query);
    }

   /**
    * Returns the values of the column `column` of the table `table`. This
    * method is equivalent to  {@link #readDoubleData(Connection,String)
    * readDoubleData(connection, "SELECT column FROM table")}.
    */
   public static double[] readDoubleData (Connection connection,
                                          String table, String column)
            throws SQLException {
        final String query = "SELECT " + column + " FROM " + table;
        
        return readDoubleData (connection, query);
    }

   /**
    * Copies the result of the SQL query `query` into an array of
    * integers. This method uses the statement `stmt` to execute the given
    * query, and assumes that the first column of the result set contains
    * integer values. Each row of the result set then becomes an element
    * of an array of integers which is returned by this method. This
    * method throws an  SQLException if the query is not valid. The given
    * statement `stmt` must not be set up to produce forward-only result
    * sets.
    *  @param stmt         the statement used to make the query.
    *  @param query        the SQL query to execute.
    *  @return the first column of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static int[] readIntData (Statement stmt, String query)
            throws SQLException {
      ResultSet rs = stmt.executeQuery (query);
      rs.last();
      int[] res = new int[rs.getRow()];
      rs.first();
        
      for (int i = 0; i < res.length; i++)
      {
         res[i] = rs.getInt (1);
         rs.next();
      }
      rs.close();
            
      return res;
   }

   /**
    * Copies the result of the SQL query `query` into an array of
    * integers. This method uses the active connection `connection` to
    * create a statement, and passes this statement, with the query, to
    * #readIntData(Statement,String), which returns an array of integers.
    *  @param connection   the active connection to the database.
    *  @param query        the SQL query to execute.
    *  @return the first column of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static int[] readIntData (Connection connection, String query)
            throws SQLException {
      Statement stmt = connection.createStatement
      (ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
      try {
         return readIntData (stmt, query);
      }
      finally {
         stmt.close();
      }
   }

   /**
    * Returns the values of the column `column` of the table `table`. This
    * method is equivalent to  {@link #readIntData(Statement,String)
    * readIntData(stmt, "SELECT column FROM table")}.
    */
   public static int[] readIntData (Statement stmt, String table,
                                    String column)
            throws SQLException {
        final String query = "SELECT " + column + " FROM " + table;
        
        return readIntData (stmt, query);
    }

   /**
    * Returns the values of the column `column` of the table `table`. This
    * method is equivalent to  {@link #readIntData(Connection,String)
    * readIntData(connection, "SELECT column FROM table")}.
    */
   public static int[] readIntData (Connection connection, String table,
                                    String column)
            throws SQLException {
        final String query = "SELECT " + column + " FROM " + table;
        
        return readIntData (connection, query);
    }

   /**
    * Copies the result of the SQL query `query` into an array of objects.
    * This method uses the statement `stmt` to execute the given query,
    * and extracts values from the first column of the obtained result set
    * by using the `getObject` method. Each row of the result set then
    * becomes an element of an array of objects which is returned by this
    * method. The type of the objects in the array depends on the column
    * type of the result set, which depends on the database and query.
    * This method throws an  SQLException if the query is not valid. The
    * given statement `stmt` must not be set up to produce forward-only
    * result sets.
    *  @param stmt         the statement used to make the query.
    *  @param query        the SQL query to execute.
    *  @return the first column of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static Object[] readObjectData (Statement stmt, String query)
            throws SQLException {
      ResultSet rs = stmt.executeQuery (query);
      rs.last();
      Object[] res = new Object[rs.getRow()];
      rs.first();
        
      for (int i = 0; i < res.length; i++)
      {
         res[i] = rs.getObject (1);
         rs.next();
      }
      rs.close();
            
      return res;
   }

   /**
    * Copies the result of the SQL query `query` into an array of objects.
    * This method uses the active connection `connection` to create a
    * statement, and passes this statement, with the query, to
    * #readObjectData(Statement,String), which returns an array of
    * integers.
    *  @param connection   the active connection to the database.
    *  @param query        the SQL query to execute.
    *  @return the first column of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static Object[] readObjectData (Connection connection,
                                          String query)
            throws SQLException {
      Statement stmt = connection.createStatement
      (ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
      try {
         return readObjectData (stmt, query);
      }
      finally {
         stmt.close();
      }
   }

   /**
    * Returns the values of the column `column` of the table `table`. This
    * method is equivalent to  {@link #readObjectData(Statement,String)
    * readObjectData(stmt, "SELECT column FROM table")}.
    */
   public static Object[] readObjectData (Statement stmt,
                                          String table, String column)
            throws SQLException {
        final String query = "SELECT " + column + " FROM " + table;
        
        return readObjectData (stmt, query);
    }

   /**
    * Returns the values of the column `column` of the table `table`. This
    * method is equivalent to  {@link #readObjectData(Connection,String)
    * readObjectData(connection, "SELECT column FROM table")}.
    */
   public static Object[] readObjectData (Connection connection,
                                          String table, String column)
            throws SQLException {
        final String query = "SELECT " + column + " FROM " + table;
        
        return readObjectData (connection, query);
    }

   /**
    * Copies the result of the SQL query `query` into a rectangular 2D
    * array of double-precision values. This method uses the statement
    * `stmt` to execute the given query, and assumes that the columns of
    * the result set contain double-precision values. Each row of the
    * result set then becomes a row of a 2D array of double-precision
    * values which is returned by this method. This method throws an
    * SQLException if the query is not valid. The given statement `stmt`
    * must not be set up to produce forward-only result sets.
    *  @param stmt         the statement used to make the query.
    *  @param query        the SQL query to execute.
    *  @return the columns of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static double[][] readDoubleData2D (Statement stmt, String query)
            throws SQLException {
      ResultSet rs = stmt.executeQuery (query);
      rs.last();
      int c = rs.getMetaData().getColumnCount();
      double[][] res = new double[rs.getRow()][c];
      rs.first();
        
      for (int i = 0; i < res.length; i++) {
         for (int j = 0; j < res[i].length; j++)
            res[i][j] = rs.getDouble (1 + j);
         rs.next();
      }
      rs.close();
            
      return res;
   }

   /**
    * Copies the result of the SQL query `query` into a rectangular 2D
    * array of double-precision values. This method uses the active
    * connection `connection` to create a statement, and passes this
    * statement, with the query, to  #readDoubleData2D(Statement,String),
    * which returns a 2D array of double-precision values.
    *  @param connection   the active connection to the database.
    *  @param query        the SQL query to execute.
    *  @return the columns of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static double[][] readDoubleData2D (Connection connection,
                                              String query)
            throws SQLException {
      Statement stmt = connection.createStatement
      (ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
      try {
         return readDoubleData2D (stmt, query);
      }
      finally {
         stmt.close();
      }
   }

   /**
    * Returns the values of the columns of the table `table`. This method
    * is equivalent to  {@link #readDoubleData2D(Statement,String)
    * readDoubleData2D(stmt, "SELECT * FROM table")}.
    */
   public static double[][] readDoubleData2DTable (Statement stmt,
                                                   String table)
            throws SQLException {
        final String query = "SELECT * FROM " + table;
        
        return readDoubleData2D (stmt, query);
    }

   /**
    * Returns the values of the columns of the table `table`. This method
    * is equivalent to  {@link #readDoubleData2D(Connection,String)
    * readDoubleData2D(connection, "SELECT * FROM table")}.
    */
   public static double[][] readDoubleData2DTable (Connection connection,
                                                   String table)
            throws SQLException {
        final String query = "SELECT * FROM " + table;
        
        return readDoubleData2D (connection, query);
    }

   /**
    * Copies the result of the SQL query `query` into a rectangular 2D
    * array of integers. This method uses the statement `stmt` to execute
    * the given query, and assumes that the columns of the result set
    * contain integers. Each row of the result set then becomes a row of a
    * 2D array of integers which is returned by this method. This method
    * throws an  SQLException if the query is not valid. The given
    * statement `stmt` must not be set up to produce forward-only result
    * sets.
    *  @param stmt         the statement used to make the query.
    *  @param query        the SQL query to execute.
    *  @return the columns of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static int[][] readIntData2D (Statement stmt, String query)
            throws SQLException {
      ResultSet rs = stmt.executeQuery (query);
      rs.last();
      int c = rs.getMetaData().getColumnCount();
      int[][] res = new int[rs.getRow()][c];
      rs.first();
        
      for (int i = 0; i < res.length; i++) {
         for (int j = 0; j < res[i].length; j++)
            res[i][j] = rs.getInt (1 + j);
         rs.next();
      }
      rs.close();
            
      return res;
   }

   /**
    * Copies the result of the SQL query `query` into a rectangular 2D
    * array of integers. This method uses the active connection
    * `connection` to create a statement, and passes this statement, with
    * the query, to  #readIntData2D(Statement,String), which returns a 2D
    * array of integers.
    *  @param connection   the active connection to the database.
    *  @param query        the SQL query to execute.
    *  @return the columns of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static int[][] readIntData2D (Connection connection, String query)
            throws SQLException {
      Statement stmt = connection.createStatement
      (ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
      try {
         return readIntData2D (stmt, query);
      }
      finally {
         stmt.close();
      }
   }

   /**
    * Returns the values of the columns of the table `table`. This method
    * is equivalent to  {@link #readIntData2D(Statement,String)
    * readIntData2D(stmt, "SELECT * FROM table")}.
    */
   public static int[][] readIntData2DTable (Statement stmt, String table)
            throws SQLException {
        final String query = "SELECT * FROM " + table;
        
        return readIntData2D (stmt, query);
    }

   /**
    * Returns the values of the columns of the table `table`. This method
    * is equivalent to  {@link #readIntData2D(Connection,String)
    * readIntData2D(connection, "SELECT * FROM table")}.
    */
   public static int[][] readIntData2DTable (Connection connection,
                                             String table)
            throws SQLException {
        final String query = "SELECT * FROM " + table;
        
        return readIntData2D (connection, query);
    }

   /**
    * Copies the result of the SQL query `query` into a rectangular 2D
    * array of objects. This method uses the statement `stmt` to execute
    * the given query, and extracts values from the obtained result set by
    * using the `getObject` method. Each row of the result set then
    * becomes a row of a 2D array of objects which is returned by this
    * method. The type of the objects in the 2D array depends on the
    * column types of the result set, which depend on the database and
    * query. This method throws an  SQLException if the query is not
    * valid. The given statement `stmt` must not be set up to produce
    * forward-only result sets.
    *  @param stmt         the statement used to make the query.
    *  @param query        the SQL query to execute.
    *  @return the columns of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static Object[][] readObjectData2D (Statement stmt, String query)
            throws SQLException {
      ResultSet rs = stmt.executeQuery (query);
      rs.last();
      int c = rs.getMetaData().getColumnCount();
      Object[][] res = new Object[rs.getRow()][c];
      rs.first();
        
      for (int i = 0; i < res.length; i++) {
         for (int j = 0; j < res[i].length; j++)
            res[i][j] = rs.getObject (1 + j);
         rs.next();
      }
      rs.close();
            
      return res;
   }

   /**
    * Copies the result of the SQL query `query` into a rectangular 2D
    * array of integers. This method uses the active connection
    * `connection` to create a statement, and passes this statement, with
    * the query, to  #readObjectData2D(Statement,String), which returns a
    * 2D array of integers.
    *  @param connection   the active connection to the database.
    *  @param query        the SQL query to execute.
    *  @return the columns of the result set.
    *
    *  @exception SQLException if the query is not valid.
    */
   public static Object[][] readObjectData2D (Connection connection,
                                              String query)
            throws SQLException {
      Statement stmt = connection.createStatement
      (ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
      try {
         return readObjectData2D (stmt, query);
      }
      finally {
         stmt.close();
      }
   }

   /**
    * Returns the values of the columns of the table `table`. This method
    * is equivalent to  {@link #readObjectData2D(Statement,String)
    * readObjectData2D(stmt, "SELECT * FROM table")}.
    */
   public static Object[][] readObjectData2DTable (Statement stmt,
                                                   String table)
            throws SQLException {
        final String query = "SELECT * FROM " + table;
        
        return readObjectData2D (stmt, query);
    }

   /**
    * Returns the values of the columns of the table `table`. This method
    * is equivalent to  {@link #readObjectData2D(Connection,String)
    * readObjectData2D(connection, "SELECT * FROM table")}.
    */
   public static Object[][] readObjectData2DTable (Connection connection,
                                                   String table)
            throws SQLException {
        final String query = "SELECT * FROM " + table;
        
        return readObjectData2D (connection, query);
    }

}